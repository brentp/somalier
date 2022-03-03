import os
import lapper
import algorithm
import hts
import slivarpkg/gnotate
import strutils
import strformat
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

var target_af = 0.48
if getEnv("SOMALIER_TARGET_AF") != "":
    target_af = parseFloat(getEnv("SOMALIER_TARGET_AF"))
    stderr.write_line &"[somalier/find-sites] using {target_af:.2f} as the target allele frequency"

proc vf_by(a:vf, b:vf): int =
  return cmp((a.af - target_af).abs, (b.af - target_af).abs)

proc start(r: region): int {.inline.} = return r.start
proc stop(r: region): int {.inline.} = return r.stop

proc bed_line_to_region(line: string): region =
  var cse = line.strip().split('\t', 5)
  if len(cse) < 3:
    stderr.write_line("[somalier] skipping bad bed line:", line.strip())
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
    let x = bed_regions.hasKeyOrPut(v.chrom, new_seq[region]())
    bed_regions[v.chrom].add(v)

  # since it is read into mem, can also well sort.
  for chrom, ivs in bed_regions.mpairs:
    sort(ivs, proc (a, b: region): int = int(a.start) - int(b.start))

  hts.free(kstr.s)

proc closest_dist(v:Variant, vs: seq[int64]): int64 =
  ## find the distance to the closest variant by brute force.
  let vstart = v.start
  result = (vstart - vs[0]).abs
  for ostart in vs:
    result = min(result, (vstart - ostart).abs)

proc findsites_main*() =
  var p = newParser("somalier find-sites"):
    option("-x", "--exclude", multiple=true, help="optional exclude files")
    option("-i", "--include", help="optional include file. only consider variants that fall in ranges within this file")
    option("--gnotate-exclude", help="sites in slivar gnotation (zip) format to exclude")
    option("--snp-dist", help="minimum distance between autosomal SNPs to avoid linkage", default="10000")
    option("--min-AN", help="minimum number of alleles (AN) at the site. (must be less than twice number of samples in the cohort)", default="115_000")
    option("--min-AF", help="minimum allele frequency for a site", default="0.15")
    arg("vcf", help="population VCF to use to find sites", nargs=1)

  var argv = commandLineParams()
  if len(argv) > 0 and argv[0] == "find-sites":
    argv = argv[1..argv.high]
  if len(argv) == 0:
    argv = @["-h"]
  let opts = p.parse(argv)
  if opts.help:
    stderr.write_line "this will write output sites to: ./sites.vcf.gz"
    quit 0

  var
    vcf:VCF
    wtr:VCF
    snp_dist = parseInt(opts.snp_dist)
    min_AN = parseInt(opts.min_AN)
    min_AF = parseFloat(opts.min_AF)
    gno:Gnotater

  var exclude_regions = newTable[string, seq[region]]()
  var include_regions = newTable[string, seq[region]]()
  var indels = newTable[string, seq[region]]()
  var snps = newTable[string, seq[region]]()
  for e in opts.exclude:
    bed_to_table(e, exclude_regions)
  if opts.include != "":
    opts.include.bed_to_table(include_regions)

  if opts.gnotate_exclude != "":
    doAssert gno.open(opts.gnotate_exclude, tmpDir=getEnv("TMPDIR"))

  if not open(vcf, opts.vcf, threads=2):
    quit "couldn't open " & opts.vcf
  vcf.set_samples(@["^"])
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
  var lap_incl:Lapper[region]
  var empty_regions: seq[region]

  var saved = newSeqOfCap[vf](100000)
  var ranksum = @[0'f32]

  for v in vcf:
    if v.c == nil: break
    if $v.CHROM notin ["chrY", "Y", "chrX", "X"] and v.REF == "C": continue
    if v.CHROM != last_chrom:
      last_chrom = v.CHROM
      stderr.write_line "on chrom:", last_chrom
      if exclude_regions.contains($last_chrom):
        lap = lapify(exclude_regions[$last_chrom])
      else:
        var seqs = newSeq[region]()
        lap = lapify(seqs)
      if include_regions.contains($last_chrom):
        lap_incl = lapify(include_regions[$last_chrom])
      else:
        var seqs = newSeq[region]()
        lap_incl = lapify(seqs)


    # stuff outside of PAR on human only.
    if $last_chrom in ["X", "chrX"] and (v.start < 2781479 or v.start > 154931044) : continue

    if v.info.get("AF", afs) != Status.OK:
      stderr.write_line "[somalier] af not found, using 0"
      afs = @[0'f32]

    if v.REF.len > 1 or v.ALT.len > 1 or v.ALT[0].len > 1:
      if afs[0] > 0.02:
        let x = indels.hasKeyOrPut($v.CHROM, newSeq[region]())
        var r = region(chrom: $v.CHROM, start: max(v.start.int - 7, 0), stop: v.stop.int + 7)
        indels[$v.CHROM].add(r)
      continue

    if v.FILTER != "PASS":
      continue

    if afs[0] > 0.01:
      var r = region(chrom: $v.CHROM, start: max(v.start.int - 1, 0), stop: v.stop.int + 1)
      #discard snps.hasKeyOrPut($v.CHROM, newSeq[region]())
      #snps[$v.CHROM].add(r)
      snps.mGetOrPut($v.CHROM, newSeq[region]()).add(r)

    # check exclude after putting into interval trees
    if gno != nil and gno.contains(v) and $v.CHROM notin ["chrX", "X"]:
      continue

    var info = v.info
    if info.get("AN", ans) == Status.OK and (($v.CHROM notin ["chrX", "X", "chrY", "Y"]) and ans[0] < min_AN) :
      continue
    if $v.CHROM in ["chrY", "Y", "chrX", "X"]:
      if afs[0] < 0.04 or afs[0] > 0.96: continue
    else:
      if afs[0] < min_AF or afs[0] > (1 - min_AF): continue

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
    if $v.CHROM notin ["chrX", "X", "chrY", "Y"]:
      var skip = false
      for rs in @["BaseQ", "MQ", "Clipping", "ReadPos"]:
        if info.get(rs & "RankSum", ranksum) == Status.OK and abs(ranksum[0]) > 2.4:
          skip = true
          break
      if skip: continue

      if info.get("FS", ranksum) == Status.OK and abs(ranksum[0]) > 2.4:
        continue
      if info.get("QD", ranksum) == Status.OK and abs(ranksum[0]) < 12:
        continue
      if info.get("MQ", ranksum) == Status.OK and ranksum[0] < 50:
        continue

    if lap.find(max(0, v.start.int - 5), v.stop.int + 5, empty_regions):
      continue

    if include_regions.len > 0 and not lap_incl.find(v.start.int, v.stop.int, empty_regions):
      continue

    #discard wtr.write_variant(v)
    var x:float32 = afs[0]
    var vfo = vf(v:v.copy(), af:x)

    saved.add(vfo)

  # wtr.close()
  if not wtr.write_header:
    quit "couldn't write vcf header"

  var indel_lappers = newTable[string, Lapper[region]]()
  for k, vs in indels.mpairs:
    indel_lappers[k] = lapify(vs)
  var snp_lappers = newTable[string, Lapper[region]]()
  for k, vs in snps.mpairs:
    snp_lappers[k] = lapify(vs)


  # Sort all variants by reverse AF, then write any variant that's
  # not within X bases of a previously added variant.
  saved.sort(vf_by)
  stderr.write_line $(saved.len), " candidate variants"
  var added = newTable[cstring, seq[int64]]()

  var used = newSeqOfCap[Variant](8192)
  var usedxy = newSeqOfCap[Variant](8192)
  var empty: seq[region]
  var xcounts = 0 # prevent too many X
  var ycounts = 0 # prevent too many y

  for v in saved:
    let chrom = $v.v.CHROM
    if chrom in indel_lappers and indel_lappers[chrom].find(max(0, v.v.start.int-1), v.v.stop.int+1, empty):
      continue
    if chrom in snp_lappers and snp_lappers[chrom].find(max(0, v.v.start.int-1), v.v.stop.int+1, empty) and empty.len > 1:
      #echo "skipping with nearby snp"
      continue
    let x = added.hasKeyOrPut(v.v.CHROM, newSeq[int64]())
    var vs = added[v.v.CHROM]

    if len(vs) > 0:
      let d = v.v.closest_dist(vs)
      if chrom in ["chrY", "Y"]:
        if d < 200:
          continue
      elif chrom in ["chrX", "X"]:
        if d < 1000:
          continue
      elif d < snp_dist:
        continue

    #discard wtr.write_variant(v.v)
    if chrom in ["chrX", "X", "Y", "chrY"]:
      if chrom in ["chrX", "X"]:
          if xcounts > 10000:
              continue
          xcounts += 1
      if chrom in ["chrY", "Y"]:
          if ycounts > 5000:
              continue
          ycounts += 1
      usedxy.add(v.v.copy())
    else:
      vs.add(v.v.start)
      used.add(v.v.copy())

    added[v.v.CHROM] = vs

  stderr.write_line &"sorted and filtered to {used.len} autosomal variants. now dropping INFOs and writing"
  if used.len.int > uint16.high.int:
    used = used[0..<uint16.high.int]

  var uu = [used, usedxy]

  for us in uu.mitems:
    sort(us, proc(a, b: Variant): int =
        if a.rid != b.rid: return cmp(a.rid, b.rid)
        return cmp(a.start, b.start)
    )
    for v in us:
      var info = v.info
      for f in info.fields:
        if f.name in ["AF", "AC"]: continue
        try:
          discard info.delete(f.name)
        except:
          discard
      # these make the compressed vcf smaller
      v.ID = "."
      v.QUAL = 100
      doAssert wtr.write_variant(v)

  stderr.write_line "[somalier] wrote ", $used.len, " variants to:", $out_path
  if vcf.header == nil:
    quit "bad"
  wtr.close()

when isMainModule:
  findsites_main()
