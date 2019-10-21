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
import ./somalierpkg/depthview
import strformat
import argparse
import strutils

proc get_sample_name(bam:Bam): string =
    var txt = newString(bam.hdr.hdr.l_text)
    copyMem(txt[0].addr, bam.hdr.hdr.text, txt.len)
    for line in txt.split("\n"):
      if line.startsWith("@RG") and "\tSM:" in line:
        result = line.split("\tSM:")[1].split("\t")[0].strip()

    if result.len == 0:
      raise newException(ValueError, "[somalier] no read-group in bam file")

proc looks_like_gvcf_variant(v:Variant): bool {.inline.} =
  result = v.c.n_allele == 1
  for a in v.ALT:
    # match either <*> or <NON_REF>, etc.
    if (a.len > 0 and a[0] == '<' and a[a.high] == '>') or a == ".":
      return true

proc get_variant(ivcf:VCF, site:Site): Variant =
  for v in ivcf.query(&"{site.chrom}:{site.position+1}-{site.position+2}"):
    if v.looks_like_gvcf_variant:
      return v.copy()

    if v.start == site.position and (
      (v.REF == site.A_allele and v.ALT[0] == site.B_allele) or
      (v.REF == site.B_allele and v.ALT[0] == site.A_allele)):
      return v.copy()

proc get_ref_alt_counts(ivcf:VCF, sites:seq[Site], fai:Fai=nil): seq[counts] =
  result = newSeq[counts](ivcf.samples.len)
  var vcf_samples = ivcf.samples
  for j, s in vcf_samples:
    result[j].sample_name = s
    result[j].sites = newSeqOfCap[allele_count](sites.len)

  var AD = newSeq[int32](5*vcf_samples.len)
  var n = 0

  for i, site in sites:
    var v = ivcf.get_variant(site)
    if v == nil or ($v.CHROM notin ["chrX", "X", "chrY", "Y"] and v.FILTER notin ["PASS", "", ".", "RefCall"]) or v.format.get("AD", AD) != Status.OK:
      AD.setLen(vcf_samples.len * 2)
      zeroMem(AD[0].addr, AD.len * sizeof(AD[0]))
      if v != nil and v.looks_like_gvcf_variant:
        var dp:seq[int32]
        if v.format.get("MIN_DP", dp) == Status.OK or v.format.get("DP", dp) == Status.OK:
           AD[0] = dp[0]
           n += 1
    else:
      n += 1
    var mult = int(AD.len / vcf_samples.len)
    for j, s in vcf_samples:
      var ac = allele_count(nref: max(0, AD[mult * j]).uint32, nalt: max(0, AD[mult * j + 1]).uint32, nother: 0)
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
    var s = newFileStream(opts.outdir & "/" & cnts.sample_name & ".somalier", fmWrite)
    s.write(formatVersion.uint8)
    s.write(cnts.sample_name.len.uint8)
    s.write(cnts.sample_name)
    s.write(cnts.sites.len.uint16)
    s.write(cnts.x_sites.len.uint16)
    s.write(cnts.y_sites.len.uint16)
    for st in cnts.sites:
      s.write(st)
    for st in cnts.x_sites:
      s.write(st)
    for st in cnts.y_sites:
      s.write(st)
    s.close()

proc main() =

  type pair = object
    f: proc()
    description: string

  var dispatcher = {
    "extract": pair(f:extract_main, description: "extract genotype-like information for a single sample from VCF/BAM/CRAM."),
    "relate": pair(f:rel_main, description: "aggregate `extract`ed information and calculate relatedness among samples."),
    "depthview": pair(f:depth_main, description: "plot per-chromosome depth for each sample for quick quality-control"),
    "find-sites": pair(f:findsites_main, description: "create a new sites.vcf.gz file from a population VCF (this is rarely needed)."),
  }.toOrderedTable


  stderr.write_line "somalier version: " & somalierVersion
  var args = commandLineParams()
  if len(args) > 0 and args[0] in dispatcher:
    dispatcher[args[0]].f()
    return

  if len(args) == 0 or args[0] in ["-h", "--help"]:
    stderr.write_line "Commands: "
    for k, v in dispatcher:
      echo &"  {k:<13}:   {v.description}"
  else:
    echo &"unknown program '{args[0]}'"
    quit ""



when isMainModule:
  main()
