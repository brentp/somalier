import os
import hts
import hts/files
import hts/private/hts_concat
import streams

import hile
import tables
import times
import ./somalierpkg/common
import ./somalierpkg/version
import ./somalierpkg/relate
import ./somalierpkg/findsites
import ./somalierpkg/ancestry
import ./somalierpkg/depthview
import strformat
import argparse
import strutils

const ENV_SAMPLE_NAME = "SOMALIER_SAMPLE_NAME"

proc get_sample_name(bam:Bam): string =
    var txt = newString(bam.hdr.hdr.l_text)
    copyMem(txt[0].addr, bam.hdr.hdr.text, txt.len)
    for line in txt.split("\n"):
      if line.startsWith("@RG") and "\tSM:" in line:
        result = line.split("\tSM:")[1].split("\t")[0].strip()

    var env_name = getEnv(ENV_SAMPLE_NAME)
    # if the result was already set, we warn.
    if result.len != 0 and env_name.len != 0:
      stderr.write_line &"[somalier] using sample name from ENV: {env_name}"
    # always use env if it was given.
    if env_name.len != 0:
      result = env_name

    if result.len == 0:
      raise newException(ValueError, &"[somalier] no read-group in bam file or given via '{ENV_SAMPLE_NAME}'")

proc looks_like_gvcf_variant(v:Variant): bool {.inline.} =
  result = v.c.n_allele == 1
  for a in v.ALT:
    # match either <*> or <NON_REF>, etc.
    if (a.len > 0 and a[0] == '<' and a[a.high] == '>') or a == ".":
      return true

proc get_variant(ivcf:VCF, site:Site): Variant =
  for v in ivcf.query(&"{site.chrom}:{site.position+1}-{site.position+2}"):
    if v.start == site.position and v.ALT.len > 0 and (
      (v.REF == site.A_allele and v.ALT[0] == site.B_allele) or
      (v.REF == site.B_allele and v.ALT[0] == site.A_allele)):
      return(v.copy())

    # we keep a gvcf variant in case we dont find an exact match.
    if v.looks_like_gvcf_variant and result == nil:
      result = v.copy()

var allowed_filters = @["PASS", "", ".", "RefCall"]
allowed_filters.add(getEnv("SOMALIER_ALLOWED_FILTERS").split(","))

proc ok(v:Variant): bool {.inline.} =
  if v == nil: return false
  if v.FILTER notin allowed_filters: return false
  return true

proc fill(AD: var seq[int32], alts: seq[int8]): bool =
  ## return true if any sample had a known genotype
  result = false
  for i, a in alts:
    case a:
      of 0:
        result = true
        AD[2 * i] = 20
        AD[2 * i + 1] = 0
      of 1:
        result = true
        AD[2 * i] = 10
        AD[2 * i + 1] = 10
      of 2:
        result = true
        AD[2 * i] = 0
        AD[2 * i + 1] = 20
      else:
        AD[2 * i] = 0
        AD[2 * i + 1] = 0


proc get_ref_alt_counts(ivcf:VCF, sites:seq[Site], fai:Fai=nil): seq[counts] =
  result = newSeq[counts](ivcf.samples.len)
  var vcf_samples = ivcf.samples
  for j, s in vcf_samples:
    result[j].sample_name = s
    result[j].sites = newSeqOfCap[allele_count](sites.len)

  var has_AD = true
  try:
    discard ivcf.header.get("AD", BCF_HEADER_TYPE.BCF_HL_FMT)
  except KeyError:
    has_AD = false

  if not has_AD:
    stderr.write_line "[somalier] FORMAT field 'AD' not found for depth information. using genotype only"

  var AD = newSeq[int32](5*vcf_samples.len)
  var x = newSeq[int32](vcf_samples.len)
  var n = 0
  var dp:seq[int32]

  for i, site in sites:
    zeroMem(AD[0].addr, AD.len * sizeof(AD[0]))
    var v = ivcf.get_variant(site)
    # NOTE that if Status is OK, then we just fill AD
    if v == nil or ($v.CHROM notin ["chrX", "X", "chrY", "Y"] and v.FILTER notin allowed_filters) or v.format.get("AD", AD) != Status.OK:
      AD.setLen(vcf_samples.len * 2)
      if v.ok and not has_AD:
        if AD.fill(v.format.genotypes(x).alts):
          n += 1
      elif v.ok and v.looks_like_gvcf_variant:
        var dp:seq[int32]
        if v.format.get("MIN_DP", dp) == Status.OK or v.format.get("DP", dp) == Status.OK:
          AD[0] = dp[0]
          n += 1
    elif v.ok and vcf_samples.len > 1 and v.looks_like_gvcf_variant:
      # we can get here with AD filled, but some samples have unknown/empty
      # values, but we can fill reference values with MIN_DP or DP depending on
      # which is available.
      var mult = v.ALT.len + 1
      if v.format.get("MIN_DP", dp) == Status.OK or v.format.get("DP", dp) == Status.OK:
        # iterate over samples and use MIN_DP/DP as the reference count
        for j in 0..<vcf_samples.len:
          if AD[j*mult] < 0:
            AD[j*mult] = max(0, dp[j])
      n += 1
    else:
      n += 1

    var mult = int(AD.len / vcf_samples.len)
    for j, s in vcf_samples:
      var ac:allele_count
      if mult >= 2:
        ac = allele_count(nref: max(0, AD[mult * j]).uint32, nalt: max(0, AD[mult * j + 1]).uint32, nother: 0)
      else:
        ac = allele_count(nref: 0'u32, nalt: 0'u32, nother: 0)

      doAssert site.A_allele < site.B_allele
      # we always store a site by alphabetical order, so sometimes we have to
      # flip the alts.
      if site.flip:
        swap(ac.nref, ac.nalt)

      case site.chrom:
        of ["X", "chrX"]:
          result[j].x_sites.add(ac)
        of ["Y", "chrY"]:
          result[j].y_sites.add(ac)
        else:
          result[j].sites.add(ac)

  stderr.write_line &"[somalier] found {n} sites"

proc get_ref_alt_counts(ibam:Bam, sites:seq[Site], fai:Fai): counts =

  result.sites = newSeqOfCap[allele_count](sites.len)
  result.sample_name = ibam.get_sample_name()

  var cfg = Config(MinMappingQuality: 1, ExcludeFlags:BAM_FUNMAP or BAM_FSECONDARY or BAM_FQCFAIL or BAM_FDUP)

  for i, site in sites:
    doAssert site.A_allele.len == 1 and site.B_allele.len == 1, "[somalier] can only genotype single nucleotode variants directly from alignments"
    var h = hileup(ibam, site.chrom, site.position, fai, cfg)
    checkSiteRef(site, fai)

    var ac = allele_count()

    for b in h.bases:
      if b.base.char == site.A_allele[0]:
        ac.nref.inc
      elif b.base.char == site.B_allele[0]:
        ac.nalt.inc
      else:
        ac.nother.inc

    case site.chrom:
      of ["X", "chrX"]:
        result.x_sites.add(ac)
      of ["Y", "chrY"]:
        result.y_sites.add(ac)
      else:
        result.sites.add(ac)


proc extract_main() =
  var argv = commandLineParams()
  if argv[0] == "extract":
    argv = argv[1..argv.high]
  if len(argv) == 0:
    argv = @["-h"]

  var p = newParser("somalier extract"):
    help("extract genotype-like information for a single-sample at selected sites")
    option("-s", "--sites", help="sites vcf file of variants to extract")
    option("-f", "--fasta", help="path to reference fasta file")
    option("-d", "--out-dir", help="path to output directory", default=".")
    option("--sample-prefix", help="prefix for the sample name stored inside the digest", default="")
    arg("sample_file", help="single-sample CRAM/BAM/GVCF file or multi/single-sample VCF from which to extract")

  let opts = p.parse(argv)
  if opts.help:
    quit 0
  if opts.sites == "":
    echo p.help
    quit "[somalier] --sites file required"
  if opts.fasta == "":
    echo p.help
    quit "[somalier] --fasta file required"

  var fai:Fai
  if not open(fai, opts.fasta):
    quit "[somalier] unable to open fasta file:" & opts.fasta

  var sites = readSites(opts.sites)
  createDir(opts.out_dir)

  var ibam: Bam
  var ivcf: VCF

  var sample_counts: seq[counts]

  if opts.sample_file.file_type in {FileType.VCF, FileType.BCF}:
    if not open(ivcf, opts.sample_file):
      quit "[somalier] couldn't open sample VCF file"
    sample_counts = ivcf.get_ref_alt_counts(sites, fai)

  else:
    if not open(ibam, opts.sample_file, fai=opts.fasta, index=true):
      quit "[somalier] couldn't open :" & opts.sample_file
    var cnts = ibam.get_ref_alt_counts(sites, fai)
    sample_counts = @[cnts]

  for cnts in sample_counts:
    var stored_sample_name : string = opts.sample_prefix & cnts.sample_name
    let fname = opts.outdir & "/" & cnts.sample_name & ".somalier"
    cnts.write_counts(stored_sample_name, fname)

proc main() =

  type pair = object
    f: proc()
    description: string

  var dispatcher = {
    "extract": pair(f:extract_main, description: "extract genotype-like information for a single sample from VCF/BAM/CRAM."),
    "relate": pair(f:rel_main, description: "aggregate `extract`ed information and calculate relatedness among samples."),
    "ancestry": pair(f:ancestry_main, description: "perform ancestry prediction on a set of samples, given a set of labeled samples"),
    #"depthview": pair(f:depth_main, description: "plot per-chromosome depth for each sample for quick quality-control"),
    "find-sites": pair(f:findsites_main, description: "create a new sites.vcf.gz file from a population VCF (this is rarely needed)."),
  }.toOrderedTable


  stderr.write_line "somalier version: " & somalierVersion
  var args = commandLineParams()
  if len(args) > 0 and args[0] in dispatcher:
    dispatcher[args[0]].f()
    return

  if len(args) == 0 or args[0] in ["-h", "--help"]:
    stdout.write_line "Commands: "
    for k, v in dispatcher:
      echo &"  {k:<13}:   {v.description}"
  else:
    echo &"unknown program '{args[0]}'"
    quit ""


when isMainModule:
  main()
