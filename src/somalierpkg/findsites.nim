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

  if argv[0] == "find-sites":
    argv = argv[1..argv.high]
  if len(argv) == 0:
    argv = @["-h"]
  let opts = p.parse(argv)
  if opts.help:
    quit 0

  var
    vcf:VCF
    wtr:VCF

  var exclude_fields = @["AS_FilterStatus", "AC","AN","DB","DP","FS","InbreedingCoeff","MQ","QD","SOR","VQSLOD","VQSR_culprit","VQSR_NEGATIVE_TRAIN_SITE","VQSR_POSITIVE_TRAIN_SITE","GQ_HIST_ALT","DP_HIST_ALT","AB_HIST_ALT","GQ_HIST_ALL","DP_HIST_ALL","AB_HIST_ALL","AC_AFR","AC_AMR","AC_ASJ","AC_EAS","AC_FIN","AC_NFE","AC_OTH","AC_SAS","AC_Male","AC_Female","AN_AFR","AN_AMR","AN_ASJ","AN_EAS","AN_FIN","AN_NFE","AN_OTH","AN_SAS","AN_Male","AN_Female","AF_AFR","AF_AMR","AF_ASJ","AF_EAS","AF_FIN","AF_NFE","AF_OTH","AF_SAS","AF_Male","AF_Female","GC_AFR","GC_AMR","GC_ASJ","GC_EAS","GC_FIN","GC_NFE","GC_OTH","GC_SAS","GC_Male","GC_Female","AC_raw","AN_raw","AF_raw","GC_raw","GC","Hom_AFR","Hom_AMR","Hom_ASJ","Hom_EAS","Hom_FIN","Hom_NFE","Hom_OTH","Hom_SAS","Hom_Male","Hom_Female","Hom_raw","Hom","STAR_AC","STAR_AC_raw","STAR_Hom","POPMAX","AC_POPMAX","AN_POPMAX","AF_POPMAX","DP_MEDIAN","DREF_MEDIAN","GQ_MEDIAN","AB_MEDIAN","AS_RF","AS_RF_POSITIVE_TRAIN","AS_RF_NEGATIVE_TRAIN","AC_AFR_Male","AC_AMR_Male","AC_ASJ_Male","AC_EAS_Male","AC_FIN_Male","AC_NFE_Male","AC_OTH_Male","AC_SAS_Male","AC_AFR_Female","AC_AMR_Female","AC_ASJ_Female","AC_EAS_Female","AC_FIN_Female","AC_NFE_Female","AC_OTH_Female","AC_SAS_Female","AN_AFR_Male","AN_AMR_Male","AN_ASJ_Male","AN_EAS_Male","AN_FIN_Male","AN_NFE_Male","AN_OTH_Male","AN_SAS_Male","AN_AFR_Female","AN_AMR_Female","AN_ASJ_Female","AN_EAS_Female","AN_FIN_Female","AN_NFE_Female","AN_OTH_Female","AN_SAS_Female","AF_AFR_Male","AF_AMR_Male","AF_ASJ_Male","AF_EAS_Male","AF_FIN_Male","AF_NFE_Male","AF_OTH_Male","AF_SAS_Male","AF_AFR_Female","AF_AMR_Female","AF_ASJ_Female","AF_EAS_Female","AF_FIN_Female","AF_NFE_Female","AF_OTH_Female","AF_SAS_Female","GC_AFR_Male","GC_AMR_Male","GC_ASJ_Male","GC_EAS_Male","GC_FIN_Male","GC_NFE_Male","GC_OTH_Male","GC_SAS_Male","GC_AFR_Female","GC_AMR_Female","GC_ASJ_Female","GC_EAS_Female","GC_FIN_Female","GC_NFE_Female","GC_OTH_Female","GC_SAS_Female","Hemi_AFR","Hemi_AMR","Hemi_ASJ","Hemi_EAS","Hemi_FIN","Hemi_NFE","Hemi_OTH","Hemi_SAS","Hemi","Hemi_raw","STAR_Hemi","CSQ","OLD_MULTIALLELIC","OLD_VARIANT","rf_tp_probability","FS","InbreedingCoeff","MQ","QD","SOR","VQSR_POSITIVE_TRAIN_SITE","VQSR_NEGATIVE_TRAIN_SITE","DP","VQSLOD","VQSR_culprit","segdup","lcr","decoy","nonpar","rf_positive_label","rf_negative_label","rf_label","rf_train","transmitted_singleton","n_alt_alleles","was_mixed","has_star","pab_max","AC_nfe_seu","AN_nfe_seu","AF_nfe_seu","nhomalt_nfe_seu","controls_AC_afr_male","controls_AN_afr_male","controls_AF_afr_male","controls_nhomalt_afr_male","non_neuro_AC_eas_kor","non_neuro_AN_eas_kor","non_neuro_AF_eas_kor","non_neuro_nhomalt_eas_kor","non_cancer_AC_asj_female","non_cancer_AN_asj_female","non_cancer_AF_asj_female","non_cancer_nhomalt_asj_female","AC_raw","AN_raw","AF_raw","nhomalt_raw","AC_fin_female","AN_fin_female","AF_fin_female","nhomalt_fin_female","non_cancer_AC_oth_female","non_cancer_AN_oth_female","non_cancer_AF_oth_female","non_cancer_nhomalt_oth_female","AC_nfe_bgr","AN_nfe_bgr","AF_nfe_bgr","nhomalt_nfe_bgr","non_neuro_AC_asj_female","non_neuro_AN_asj_female","non_neuro_AF_asj_female","non_neuro_nhomalt_asj_female","AC_sas_male","AN_sas_male","AF_sas_male","nhomalt_sas_male","non_neuro_AC_afr_male","non_neuro_AN_afr_male","non_neuro_AF_afr_male","non_neuro_nhomalt_afr_male","AC_afr_male","AN_afr_male","AF_afr_male","nhomalt_afr_male","AC_afr","AN_afr","AF_afr","nhomalt_afr","controls_AC_nfe_swe","controls_AN_nfe_swe","controls_AF_nfe_swe","controls_nhomalt_nfe_swe","non_neuro_AC_afr_female","non_neuro_AN_afr_female","non_neuro_AF_afr_female","non_neuro_nhomalt_afr_female","non_cancer_AC_female","non_cancer_AN_female","non_cancer_AF_female","non_cancer_nhomalt_female","non_cancer_AC_nfe_onf","non_cancer_AN_nfe_onf","non_cancer_AF_nfe_onf","non_cancer_nhomalt_nfe_onf","non_cancer_AC_male","non_cancer_AN_male","non_cancer_AF_male","non_cancer_nhomalt_male","AC_eas_female","AN_eas_female","AF_eas_female","nhomalt_eas_female","non_cancer_AC_sas_female","non_cancer_AN_sas_female","non_cancer_AF_sas_female","non_cancer_nhomalt_sas_female","AC_afr_female","AN_afr_female","AF_afr_female","nhomalt_afr_female","AC_sas","AN_sas","AF_sas","nhomalt_sas","non_neuro_AC_female","non_neuro_AN_female","non_neuro_AF_female","non_neuro_nhomalt_female","controls_AC_afr","controls_AN_afr","controls_AF_afr","controls_nhomalt_afr","non_neuro_AC_eas_jpn","non_neuro_AN_eas_jpn","non_neuro_AF_eas_jpn","non_neuro_nhomalt_eas_jpn","AC_nfe_onf","AN_nfe_onf","AF_nfe_onf","nhomalt_nfe_onf","non_cancer_AC_amr_male","non_cancer_AN_amr_male","non_cancer_AF_amr_male","non_cancer_nhomalt_amr_male","controls_AC_fin_male",
    "controls_AN_fin_male","controls_AF_fin_male","controls_nhomalt_fin_male","non_neuro_AC_nfe_nwe","non_neuro_AN_nfe_nwe","non_neuro_AF_nfe_nwe","non_neuro_nhomalt_nfe_nwe","AC_fin_male","AN_fin_male","AF_fin_male","nhomalt_fin_male","AC_nfe_female","AN_nfe_female","AF_nfe_female","nhomalt_nfe_female","AC_amr","AN_amr","AF_amr","nhomalt_amr","non_neuro_AC_sas","non_neuro_AN_sas","non_neuro_AF_sas","non_neuro_nhomalt_sas","non_cancer_AC_fin_male","non_cancer_AN_fin_male","non_cancer_AF_fin_male","non_cancer_nhomalt_fin_male","non_cancer_AC_nfe_seu","non_cancer_AN_nfe_seu","non_cancer_AF_nfe_seu","non_cancer_nhomalt_nfe_seu","AC_eas","AN_eas","AF_eas","nhomalt_eas","nhomalt","non_neuro_AC_nfe_female","non_neuro_AN_nfe_female","non_neuro_AF_nfe_female","non_neuro_nhomalt_nfe_female","non_neuro_AC_afr","non_neuro_AN_afr","non_neuro_AF_afr","non_neuro_nhomalt_afr","controls_AC_raw","controls_AN_raw","controls_AF_raw","controls_nhomalt_raw","non_cancer_AC_eas","non_cancer_AN_eas","non_cancer_AF_eas","non_cancer_nhomalt_eas","non_cancer_AC_amr_female","non_cancer_AN_amr_female","non_cancer_AF_amr_female","non_cancer_nhomalt_amr_female","non_neuro_AC_nfe_swe","non_neuro_AN_nfe_swe","non_neuro_AF_nfe_swe","non_neuro_nhomalt_nfe_swe","controls_AC_male","controls_AN_male","controls_AF_male","controls_nhomalt_male","controls_AC_eas_jpn","controls_AN_eas_jpn","controls_AF_eas_jpn","controls_nhomalt_eas_jpn","controls_AC_nfe_female","controls_AN_nfe_female","controls_AF_nfe_female","controls_nhomalt_nfe_female","non_neuro_AC_amr","non_neuro_AN_amr","non_neuro_AF_amr","non_neuro_nhomalt_amr","non_neuro_AC_eas_female","non_neuro_AN_eas_female","non_neuro_AF_eas_female","non_neuro_nhomalt_eas_female","AC_asj_male","AN_asj_male","AF_asj_male","nhomalt_asj_male","controls_AC_nfe_male","controls_AN_nfe_male","controls_AF_nfe_male","controls_nhomalt_nfe_male","non_neuro_AC_fin","non_neuro_AN_fin","non_neuro_AF_fin","non_neuro_nhomalt_fin","non_cancer_AC_nfe_female","non_cancer_AN_nfe_female","non_cancer_AF_nfe_female","non_cancer_nhomalt_nfe_female","AC_oth_female","AN_oth_female","AF_oth_female","nhomalt_oth_female","non_cancer_AC_asj","non_cancer_AN_asj","non_cancer_AF_asj","non_cancer_nhomalt_asj","AC_nfe_swe","AN_nfe_swe","AF_nfe_swe","nhomalt_nfe_swe","controls_AC_nfe","controls_AN_nfe","controls_AF_nfe","controls_nhomalt_nfe","controls_AC_oth_female","controls_AN_oth_female","controls_AF_oth_female","controls_nhomalt_oth_female","controls_AC_asj","controls_AN_asj","controls_AF_asj","controls_nhomalt_asj","non_neuro_AC_amr_male","non_neuro_AN_amr_male","non_neuro_AF_amr_male","non_neuro_nhomalt_amr_male","controls_AC_nfe_nwe","controls_AN_nfe_nwe","controls_AF_nfe_nwe","controls_nhomalt_nfe_nwe","AC_nfe_nwe","AN_nfe_nwe","AF_nfe_nwe","nhomalt_nfe_nwe","controls_AC_nfe_seu","controls_AN_nfe_seu","controls_AF_nfe_seu","controls_nhomalt_nfe_seu","controls_AC_sas_female","controls_AN_sas_female","controls_AF_sas_female","controls_nhomalt_sas_female","non_neuro_AC_amr_female","non_neuro_AN_amr_female","non_neuro_AF_amr_female","non_neuro_nhomalt_amr_female","non_cancer_AC_eas_jpn","non_cancer_AN_eas_jpn","non_cancer_AF_eas_jpn","non_cancer_nhomalt_eas_jpn","non_neuro_AC_nfe_onf","non_neuro_AN_nfe_onf","non_neuro_AF_nfe_onf","non_neuro_nhomalt_nfe_onf","AC_eas_jpn","AN_eas_jpn","AF_eas_jpn","nhomalt_eas_jpn","non_cancer_AC_afr_male","non_cancer_AN_afr_male","non_cancer_AF_afr_male","non_cancer_nhomalt_afr_male","non_cancer_AC_afr","non_cancer_AN_afr","non_cancer_AF_afr","non_cancer_nhomalt_afr","controls_AC_amr_female","controls_AN_amr_female","controls_AF_amr_female","controls_nhomalt_amr_female","non_neuro_AC_fin_male","non_neuro_AN_fin_male","non_neuro_AF_fin_male","non_neuro_nhomalt_fin_male","AC_female","AN_female","AF_female","nhomalt_female","non_neuro_AC_nfe_bgr","non_neuro_AN_nfe_bgr","non_neuro_AF_nfe_bgr","non_neuro_nhomalt_nfe_bgr","non_neuro_AC_oth_male","non_neuro_AN_oth_male","non_neuro_AF_oth_male","non_neuro_nhomalt_oth_male","non_cancer_AC_amr","non_cancer_AN_amr","non_cancer_AF_amr","non_cancer_nhomalt_amr","controls_AC_eas_kor","controls_AN_eas_kor","controls_AF_eas_kor","controls_nhomalt_eas_kor","controls_AC_eas_male","controls_AN_eas_male","controls_AF_eas_male","controls_nhomalt_eas_male","controls_AC_oth_male","controls_AN_oth_male","controls_AF_oth_male","controls_nhomalt_oth_male","controls_AC_fin","controls_AN_fin","controls_AF_fin","controls_nhomalt_fin","AC_eas_kor","AN_eas_kor","AF_eas_kor","nhomalt_eas_kor","non_neuro_AC_nfe","non_neuro_AN_nfe","non_neuro_AF_nfe","non_neuro_nhomalt_nfe","non_neuro_AC_fin_female","non_neuro_AN_fin_female","non_neuro_AF_fin_female","non_neuro_nhomalt_fin_female","non_cancer_AC_nfe_male","non_cancer_AN_nfe_male","non_cancer_AF_nfe_male","non_cancer_nhomalt_nfe_male","controls_AC_eas_oea","controls_AN_eas_oea","controls_AF_eas_oea","controls_nhomalt_eas_oea","controls_AC_eas_female","controls_AN_eas_female","controls_AF_eas_female","controls_nhomalt_eas_female","controls_AC_nfe_onf","controls_AN_nfe_onf","controls_AF_nfe_onf","controls_nhomalt_nfe_onf","non_neuro_AC","non_neuro_AN","non_neuro_AF","non_neuro_nhomalt","AC_eas_oea","AN_eas_oea","AF_eas_oea","nhomalt_eas_oea","non_cancer_AC_oth","non_cancer_AN_oth","non_cancer_AF_oth","non_cancer_nhomalt_oth","non_neuro_AC_nfe_est","non_neuro_AN_nfe_est","non_neuro_AF_nfe_est","non_neuro_nhomalt_nfe_est","non_cancer_AC_oth_male","non_cancer_AN_oth_male","non_cancer_AF_oth_male","non_cancer_nhomalt_oth_male","AC_nfe_est","AN_nfe_est","AF_nfe_est","nhomalt_nfe_est","non_cancer_AC_afr_female","non_cancer_AN_afr_female","non_cancer_AF_afr_female","non_cancer_nhomalt_afr_female","AC_eas_male","AN_eas_male","AF_eas_male","nhomalt_eas_male","controls_AC_eas","controls_AN_eas","controls_AF_eas","controls_nhomalt_eas","non_neuro_AC_eas_male","non_neuro_AN_eas_male","non_neuro_AF_eas_male","non_neuro_nhomalt_eas_male","non_cancer_AC_nfe_nwe","non_cancer_AN_nfe_nwe","non_cancer_AF_nfe_nwe","non_cancer_nhomalt_nfe_nwe","controls_AC_sas","controls_AN_sas","controls_AF_sas","controls_nhomalt_sas","non_neuro_AC_sas_male","non_neuro_AN_sas_male","non_neuro_AF_sas_male","non_neuro_nhomalt_sas_male","non_neuro_AC_asj_male","non_neuro_AN_asj_male","non_neuro_AF_asj_male","non_neuro_nhomalt_asj_male","non_cancer_AC_nfe_bgr","non_cancer_AN_nfe_bgr","non_cancer_AF_nfe_bgr","non_cancer_nhomalt_nfe_bgr","controls_AC_oth","controls_AN_oth","controls_AF_oth","controls_nhomalt_oth","non_cancer_AC_eas_female","non_cancer_AN_eas_female","non_cancer_AF_eas_female","non_cancer_nhomalt_eas_female","AC_nfe","AN_nfe","AF_nfe","nhomalt_nfe","non_neuro_AC_asj","non_neuro_AN_asj","non_neuro_AF_asj","non_neuro_nhomalt_asj","non_neuro_AC_raw",
    "non_neuro_AN_raw","non_neuro_AF_raw","non_neuro_nhomalt_raw","non_cancer_AC_asj_male","non_cancer_AN_asj_male","non_cancer_AF_asj_male","non_cancer_nhomalt_asj_male","AC_fin","AN_fin","AF_fin","nhomalt_fin","AC_nfe_male","AN_nfe_male","AF_nfe_male","nhomalt_nfe_male","controls_AC_amr_male","controls_AN_amr_male","controls_AF_amr_male","controls_nhomalt_amr_male","non_neuro_AC_eas_oea","non_neuro_AN_eas_oea","non_neuro_AF_eas_oea","non_neuro_nhomalt_eas_oea","AC_sas_female","AN_sas_female","AF_sas_female","nhomalt_sas_female","controls_AC_afr_female","controls_AN_afr_female","controls_AF_afr_female","controls_nhomalt_afr_female","controls_AC_amr","controls_AN_amr","controls_AF_amr","controls_nhomalt_amr","AC_asj_female","AN_asj_female","AF_asj_female","nhomalt_asj_female","non_cancer_AC_nfe_est","non_cancer_AN_nfe_est","non_cancer_AF_nfe_est","non_cancer_nhomalt_nfe_est","non_neuro_AC_eas","non_neuro_AN_eas","non_neuro_AF_eas","non_neuro_nhomalt_eas","non_cancer_AC_nfe","non_cancer_AN_nfe","non_cancer_AF_nfe","non_cancer_nhomalt_nfe","non_neuro_AC_male","non_neuro_AN_male","non_neuro_AF_male","non_neuro_nhomalt_male","non_neuro_AC_sas_female","non_neuro_AN_sas_female","non_neuro_AF_sas_female","non_neuro_nhomalt_sas_female","AC_asj","AN_asj","AF_asj","nhomalt_asj","controls_AC_nfe_est","controls_AN_nfe_est","controls_AF_nfe_est","controls_nhomalt_nfe_est","non_cancer_AC_nfe_swe","non_cancer_AN_nfe_swe","non_cancer_AF_nfe_swe","non_cancer_nhomalt_nfe_swe","non_cancer_AC","non_cancer_AN","non_cancer_AF","non_cancer_nhomalt","non_cancer_AC_fin_female","non_cancer_AN_fin_female","non_cancer_AF_fin_female","non_cancer_nhomalt_fin_female","AC_oth","AN_oth","AF_oth","nhomalt_oth","non_neuro_AC_nfe_male","non_neuro_AN_nfe_male","non_neuro_AF_nfe_male","non_neuro_nhomalt_nfe_male","controls_AC_female","controls_AN_female","controls_AF_female","controls_nhomalt_female","non_cancer_AC_fin","non_cancer_AN_fin","non_cancer_AF_fin","non_cancer_nhomalt_fin","non_cancer_AC_eas_oea","non_cancer_AN_eas_oea","non_cancer_AF_eas_oea","non_cancer_nhomalt_eas_oea","non_cancer_AC_sas_male","non_cancer_AN_sas_male","non_cancer_AF_sas_male","non_cancer_nhomalt_sas_male","controls_AC_asj_male","controls_AN_asj_male","controls_AF_asj_male","controls_nhomalt_asj_male","non_cancer_AC_raw","non_cancer_AN_raw","non_cancer_AF_raw","non_cancer_nhomalt_raw","non_cancer_AC_eas_male","non_cancer_AN_eas_male","non_cancer_AF_eas_male","non_cancer_nhomalt_eas_male","non_neuro_AC_oth","non_neuro_AN_oth","non_neuro_AF_oth","non_neuro_nhomalt_oth","AC_male","AN_male","AF_male","nhomalt_male","controls_AC_fin_female","controls_AN_fin_female","controls_AF_fin_female","controls_nhomalt_fin_female","controls_AC_nfe_bgr","controls_AN_nfe_bgr","controls_AF_nfe_bgr","controls_nhomalt_nfe_bgr","controls_AC_asj_female","controls_AN_asj_female","controls_AF_asj_female","controls_nhomalt_asj_female","AC_amr_male","AN_amr_male","AF_amr_male","nhomalt_amr_male","AC_amr_female","AN_amr_female","AF_amr_female","nhomalt_amr_female","AC_oth_male","AN_oth_male","AF_oth_male","nhomalt_oth_male","non_cancer_AC_sas","non_cancer_AN_sas","non_cancer_AF_sas","non_cancer_nhomalt_sas","non_neuro_AC_nfe_seu","non_neuro_AN_nfe_seu","non_neuro_AF_nfe_seu","non_neuro_nhomalt_nfe_seu","non_cancer_AC_eas_kor","non_cancer_AN_eas_kor","non_cancer_AF_eas_kor","non_cancer_nhomalt_eas_kor","controls_AC_sas_male","controls_AN_sas_male","controls_AF_sas_male","controls_nhomalt_sas_male","controls_AC","controls_AN","controls_AF","controls_nhomalt","non_neuro_AC_oth_female","non_neuro_AN_oth_female","non_neuro_AF_oth_female","non_neuro_nhomalt_oth_female","popmax","AC_popmax","AN_popmax","AF_popmax","nhomalt_popmax","non_neuro_popmax","non_neuro_AC_popmax","non_neuro_AN_popmax","non_neuro_AF_popmax","non_neuro_nhomalt_popmax","non_cancer_popmax","non_cancer_AC_popmax","non_cancer_AN_popmax","non_cancer_AF_popmax","non_cancer_nhomalt_popmax","controls_popmax","controls_AC_popmax","controls_AN_popmax","controls_AF_popmax","controls_nhomalt_popmax"]


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
    if $v.CHROM in ["chrY", "Y"]:
      if afs[0] < 0.04 or afs[0] > 0.96: continue
    else:
      if afs[0] < 0.08 or afs[0] > 0.92: continue

    if info.get("AS_FilterStatus", oms) == Status.OK and oms != "PASS":
      continue
    if info.get("OLD_MULTIALLELIC", oms) == Status.OK:
      continue
    if info.get("OLD_VARIANT", oms) == Status.OK:
      continue
    if info.has_flag("segdup"):
      continue
    if info.has_flag("lcr"):
      continue

    # BaseQRankSum=0.571;ClippingRankSum=0;MQRankSum=0.101;ReadPosRankSum
    var skip = false
    for rs in @["BaseQ", "MQ", "Clipping", "ReadPos"]:
      if info.get(rs & "RankSum", ranksum) == Status.OK and abs(ranksum[0]) > 3:
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
      elif (close - v.v.start).abs < 12000:
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
    for f in exclude_fields:
      try:
        discard info.delete(f)
      except:
        discard
    doAssert wtr.write_variant(v)

  stderr.write_line "[somalier] wrote ", $used.len, " variants to:", $out_path
  if vcf.header == nil:
    quit "bad"
  wtr.close()

when isMainModule:
  findsites_main()
