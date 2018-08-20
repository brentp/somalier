import os
import lapper
import algorithm
import hts
import strutils
import tables


type region_t = ref object
  chrom: string
  start: int
  stop: int
  name: string

type vf = object
  v: Variant
  af: float32
  mqrs: float32

const target_af = 0.5

proc vf_by(a:vf, b:vf): int =
  ## sort close to the given AF.
  if (a.af - b.af) < 0.01:
    return cmp(a.mqrs.abs, b.mqrs.abs)
  return cmp((a.af - target_af).abs, (b.af - target_af).abs)

proc start(r: region_t): int {.inline.} = return r.start
proc stop(r: region_t): int {.inline.} = return r.stop


proc bed_line_to_region(line: string): region_t =
  var cse = line.strip().split('\t', 5)
  if len(cse) < 3:
    stderr.write_line("[mosdepth] skipping bad bed line:", line.strip())
    return nil
  var
    s = parse_int(cse[1])
    e = parse_int(cse[2])
    reg = region_t(chrom: cse[0], start: s, stop: e)
  if len(cse) > 3:
    reg.name = cse[3]
  return reg


proc bed_to_table(bed: string, bed_regions:var TableRef[string, seq[region_t]]) =
  var kstr = kstring_t(l:0, m: 0, s: nil)
  var hf = hts_open(cstring(bed), "r")
  while hts_getline(hf, cint(10), addr kstr) > 0:
    if ($kstr.s).startswith("track "):
      continue
    if $kstr.s[0] == "#":
      continue
    var v = bed_line_to_region($kstr.s)
    if v == nil: continue
    discard bed_regions.hasKeyOrPut(v.chrom, new_seq[region_t]())
    bed_regions[v.chrom].add(v)

  # since it is read into mem, can also well sort.
  for chrom, ivs in bed_regions.mpairs:
    sort(ivs, proc (a, b: region_t): int = int(a.start) - int(b.start))

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

when isMainModule:
  var 
    params = commandLineParams()
    vcf_path = params[0]
    vcf:VCF
    wtr:VCF


  var exclude_fields = @["AS_FilterStatus", "AC","AN","DB","DP","FS","InbreedingCoeff","MQ","QD","SOR","VQSLOD","VQSR_culprit","VQSR_NEGATIVE_TRAIN_SITE","VQSR_POSITIVE_TRAIN_SITE","GQ_HIST_ALT","DP_HIST_ALT","AB_HIST_ALT","GQ_HIST_ALL","DP_HIST_ALL","AB_HIST_ALL","AC_AFR","AC_AMR","AC_ASJ","AC_EAS","AC_FIN","AC_NFE","AC_OTH","AC_SAS","AC_Male","AC_Female","AN_AFR","AN_AMR","AN_ASJ","AN_EAS","AN_FIN","AN_NFE","AN_OTH","AN_SAS","AN_Male","AN_Female","AF_AFR","AF_AMR","AF_ASJ","AF_EAS","AF_FIN","AF_NFE","AF_OTH","AF_SAS","AF_Male","AF_Female","GC_AFR","GC_AMR","GC_ASJ","GC_EAS","GC_FIN","GC_NFE","GC_OTH","GC_SAS","GC_Male","GC_Female","AC_raw","AN_raw","AF_raw","GC_raw","GC","Hom_AFR","Hom_AMR","Hom_ASJ","Hom_EAS","Hom_FIN","Hom_NFE","Hom_OTH","Hom_SAS","Hom_Male","Hom_Female","Hom_raw","Hom","STAR_AC","STAR_AC_raw","STAR_Hom","POPMAX","AC_POPMAX","AN_POPMAX","AF_POPMAX","DP_MEDIAN","DREF_MEDIAN","GQ_MEDIAN","AB_MEDIAN","AS_RF","AS_RF_POSITIVE_TRAIN","AS_RF_NEGATIVE_TRAIN","AC_AFR_Male","AC_AMR_Male","AC_ASJ_Male","AC_EAS_Male","AC_FIN_Male","AC_NFE_Male","AC_OTH_Male","AC_SAS_Male","AC_AFR_Female","AC_AMR_Female","AC_ASJ_Female","AC_EAS_Female","AC_FIN_Female","AC_NFE_Female","AC_OTH_Female","AC_SAS_Female","AN_AFR_Male","AN_AMR_Male","AN_ASJ_Male","AN_EAS_Male","AN_FIN_Male","AN_NFE_Male","AN_OTH_Male","AN_SAS_Male","AN_AFR_Female","AN_AMR_Female","AN_ASJ_Female","AN_EAS_Female","AN_FIN_Female","AN_NFE_Female","AN_OTH_Female","AN_SAS_Female","AF_AFR_Male","AF_AMR_Male","AF_ASJ_Male","AF_EAS_Male","AF_FIN_Male","AF_NFE_Male","AF_OTH_Male","AF_SAS_Male","AF_AFR_Female","AF_AMR_Female","AF_ASJ_Female","AF_EAS_Female","AF_FIN_Female","AF_NFE_Female","AF_OTH_Female","AF_SAS_Female","GC_AFR_Male","GC_AMR_Male","GC_ASJ_Male","GC_EAS_Male","GC_FIN_Male","GC_NFE_Male","GC_OTH_Male","GC_SAS_Male","GC_AFR_Female","GC_AMR_Female","GC_ASJ_Female","GC_EAS_Female","GC_FIN_Female","GC_NFE_Female","GC_OTH_Female","GC_SAS_Female","Hemi_AFR","Hemi_AMR","Hemi_ASJ","Hemi_EAS","Hemi_FIN","Hemi_NFE","Hemi_OTH","Hemi_SAS","Hemi","Hemi_raw","STAR_Hemi","CSQ","OLD_MULTIALLELIC","OLD_VARIANT"]


  var exclude_regions = newTable[string, seq[region_t]]()
  for i in 1..<params.len:
    bed_to_table(params[i], exclude_regions)

  if not open(vcf, vcf_path, threads=2):
    quit "couldn't open " & vcf_path
  if not open(wtr, "-", mode="wz", threads=1):
    quit "couldn't open stdout for writing sites"
  wtr.header = vcf.header

  var
    last_chrom = "".cstring
    last_pos = -2000
  var afs = newSeq[float32](1)
  var oms = ""
  var lap:Lapper[region_t]
  var empty_regions: seq[region_t]

  var saved = newSeqOfCap[vf](100000)

  for v in vcf:
    if v.REF == "C": continue
    if v.CHROM != last_chrom:
      last_chrom = v.CHROM
      last_pos = -2000
      stderr.write_line $last_chrom
      if exclude_regions.contains($last_chrom):
        lap = lapify(exclude_regions[$last_chrom])
      else:
        var seqs = newSeq[region_t]()
        lap = lapify(seqs)
    if v.start - last_pos < 2000:
      continue
    if $last_chrom == "X" or $last_chrom == "Y": continue
    if $last_chrom == "chrX" or $last_chrom == "chrY": continue

    if v.REF.len > 1 or v.ALT.len > 1 or v.ALT[0].len > 1: continue
    var f = v.FILTER
    #if f != "PASS" and not ("RF" in f or "AC0" in f or "LCR" in f or "SEGDUP" in f or "InbreedingCoeff" in f):
    #  quit v.FILTER
    if f != "PASS":
      continue

    var info = v.info
    if info.floats("AF", afs) != Status.OK:
      continue
    if info.strings("AS_FilterStatus", oms) == Status.OK and oms != "PASS":
      continue
    if info.strings("OLD_MULTIALLELIC", oms) == Status.OK:
      continue
    if info.strings("OLD_VARIANT", oms) == Status.OK:
      continue
    if afs[0] < 0.08 or afs[0] > 0.92: continue

    if lap.find(v.start, v.stop, empty_regions):
      continue

    for f in exclude_fields:
      try:
        discard info.delete(f)
      except:
        discard
    #discard wtr.write_variant(v)
    #last_pos = v.start
    var x:float32 = afs[0]
    var vfo = vf(v:v.copy(), af:x)
    if info.floats("MQRankSum", afs) == Status.OK:
      vfo.mqrs = afs[0]
    else:
      stderr.write_line "no MQRankSum"

    saved.add(vfo)

  # wtr.close()
  if not wtr.write_header:
    quit "couldn't write vcf header"

  # Sort all variants by reverse AF, then write any variant that's 
  # not within X bases of a previously added variant.
  saved.sort(vf_by)
  stderr.write_line $(saved.len), " candidate variants"
  var added = newTable[cstring, seq[Variant]]()
  var wrote = 0

  for v in saved:
    discard added.hasKeyOrPut(v.v.CHROM, newSeq[Variant]())
    var vs = added.mget(v.v.CHROM)

    if len(vs) > 0 and (v.v.closest(vs).start - v.v.start).abs < 1000:
      continue

    discard wtr.write_variant(v.v)
    vs.add(v.v.copy())
    wrote += 1

    added[v.v.CHROM] = vs

  stderr.write_line "wrote ", $(wrote), " variants"
  if vcf.header == nil:
    quit "bad"
  wtr.close()
