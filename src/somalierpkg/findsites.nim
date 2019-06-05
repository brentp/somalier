import os
import lapper
import algorithm
import hts
import strutils
import tables
import argparse


type region = object
  chrom: string
  start: int
  stop: int
  name: string

type vf = object
  v: Variant
  af: float32

const target_af = 0.48

proc vf_by(a:vf, b:vf): int =
  return cmp((a.af - target_af).abs, (b.af - target_af).abs)

proc start(r: region): int {.inline.} = return r.start
proc stop(r: region): int {.inline.} = return r.stop

proc bed_line_to_region(line: string): region =
  var cse = line.strip().split('\t', 5)
  if len(cse) < 3:
    stderr.write_line("[mosdepth] skipping bad bed line:", line.strip())
    return
  var
    s = parse_int(cse[1])
    e = parse_int(cse[2])
    reg = region(chrom: cse[0], start: s, stop: e)
  if len(cse) > 3:
    reg.name = cse[3]
  return reg


proc bed_to_table(bed: string, bed_regions:var TableRef[string, seq[region]]) =
  var kstr = kstring_t(l:0, m: 0, s: nil)
  var hf = hts_open(cstring(bed), "r")
  while hts_getline(hf, cint(10), addr kstr) > 0:
    if ($kstr.s).startswith("track "):
      continue
    if $kstr.s[0] == "#":
      continue
    var v = bed_line_to_region($kstr.s)
    if v.chrom.len == 0: continue
    discard bed_regions.hasKeyOrPut(v.chrom, new_seq[region]())
    bed_regions[v.chrom].add(v)

  # since it is read into mem, can also well sort.
  for chrom, ivs in bed_regions.mpairs:
    sort(ivs, proc (a, b: region): int = int(a.start) - int(b.start))

  hts.free(kstr.s)

proc closest(v:Variant, vs: seq[Variant]): Variant =
  ## find the closest variant by brute force.
  result = vs[0]
  var closest_dist = (v.start - vs[0].start).abs
  for o in vs:
    var d = (v.start - o.start).abs
    if d < closest_dist:
      closest_dist = d
      result = o
proc findsites_main*() =
  var p = newParser("somalier find-sites"):
    option("-x", "--exclude", multiple=true, help="optional exclude files")
    arg("vcf", help="population VCF to use to find sites", nargs=1)

  var argv = commandLineParams()
  if len(argv) == 0:
    argv = @["-h"]
  elif argv[0] == "find-sites":
    argv = argv[1..argv.high]
  let opts = p.parse(argv)
  if opts.help:
    quit 0

  var
    vcf:VCF
    wtr:VCF

  var exclude_regions = newTable[string, seq[region]]()
  for e in opts.exclude:
    bed_to_table(e, exclude_regions)

  if not open(vcf, opts.vcf, threads=2):
    quit "couldn't open " & opts.vcf
  var out_path = "sites.vcf.gz"
  if not open(wtr, out_path, mode="wz", threads=1):
    quit "couldn't open stdout for writing sites"
  wtr.header = vcf.header

  var
    last_chrom = "".cstring
  var afs = newSeq[float32](1)
  var ans = newSeq[int32](1)
  var oms = ""
  var lap:Lapper[region]
  var empty_regions: seq[region]

  var saved = newSeqOfCap[vf](100000)
  var ranksum = @[0'f32]

  for v in vcf:
    if $v.CHROM notin ["chrY", "Y", "chrX", "X"] and v.REF == "C": continue
    if v.CHROM != last_chrom:
      last_chrom = v.CHROM
      if exclude_regions.contains($last_chrom):
        lap = lapify(exclude_regions[$last_chrom])
      else:
        var seqs = newSeq[region]()
        lap = lapify(seqs)

    # stuff outside of PAR on human only.
    if $last_chrom in ["X", "chrX"] and (v.start < 2781479 or v.start > 154931044) : continue

    if v.REF.len > 1 or v.ALT.len > 1 or v.ALT[0].len > 1: continue
    var f = v.FILTER
    #if f != "PASS" and not ("RF" in f or "AC0" in f or "LCR" in f or "SEGDUP" in f or "InbreedingCoeff" in f):
    #  quit v.FILTER
    if f != "PASS":
      continue

    var info = v.info
    if info.get("AF", afs) != Status.OK:
      continue
    if info.get("AN", ans) == Status.OK and (($v.CHROM notin ["chrX", "X", "chrY", "Y"]) and ans[0] < 115_000) :
      continue
    if $v.CHROM in ["chrY", "Y"]:
      if afs[0] < 0.04 or afs[0] > 0.96: continue
    else:
      if afs[0] < 0.12 or afs[0] > 0.88: continue

    if info.get("AS_FilterStatus", oms) == Status.OK and oms != "PASS":
      continue
    if info.get("OLD_MULTIALLELIC", oms) == Status.OK:
      continue
    if info.get("OLD_VARIANT", oms) == Status.OK:
      continue

    if $v.CHROM notin ["chrY", "Y"] and info.has_flag("segdup"):
      continue
    if info.has_flag("lcr"):
      continue

    # BaseQRankSum=0.571;ClippingRankSum=0;MQRankSum=0.101;ReadPosRankSum
    var skip = false
    for rs in @["BaseQ", "MQ", "Clipping", "ReadPos"]:
      if info.get(rs & "RankSum", ranksum) == Status.OK and abs(ranksum[0]) > 2.8:
        skip = true
        break
    if skip: continue

    if lap.find(v.start, v.stop, empty_regions):
      continue

    #discard wtr.write_variant(v)
    var x:float32 = afs[0]
    var vfo = vf(v:v.copy(), af:x)

    saved.add(vfo)

  # wtr.close()
  if not wtr.write_header:
    quit "couldn't write vcf header"

  # Sort all variants by reverse AF, then write any variant that's
  # not within X bases of a previously added variant.
  saved.sort(vf_by)
  stderr.write_line $(saved.len), " candidate variants"
  var added = newTable[cstring, seq[Variant]]()

  var used = newSeqOfCap[Variant](8192)

  for v in saved:
    discard added.hasKeyOrPut(v.v.CHROM, newSeq[Variant]())
    var vs = added[v.v.CHROM]

    if len(vs) > 0:
      var close = v.v.closest(vs).start
      if $v.v.CHROM in ["chrY", "Y"]:
        if (close - v.v.start).abs < 200:
          continue
      elif $v.v.CHROM in ["chrX", "X"]:
        if (close - v.v.start).abs < 20000:
          continue
      elif (close - v.v.start).abs < 10_000:
        continue

    #discard wtr.write_variant(v.v)
    vs.add(v.v.copy())
    used.add(v.v.copy())

    added[v.v.CHROM] = vs

  sort(used, proc(a, b: Variant): int =
      if a.rid != b.rid: return cmp(a.rid, b.rid)
      return cmp(a.start, b.start)
  )
  for v in used:
    var info = v.info
    for f in info.fields:
      if f.name in ["AF", "AC"]: continue
      try:
        discard info.delete(f.name)
      except:
        discard
    doAssert wtr.write_variant(v)

  stderr.write_line "[somalier] wrote ", $used.len, " variants to:", $out_path
  if vcf.header == nil:
    quit "bad"
  wtr.close()

when isMainModule:
  findsites_main()
