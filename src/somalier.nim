import os
import hts
import hts/private/hts_concat
import streams

import hile
import tables
import times
import ./somalierpkg/version
import ./somalierpkg/relate
import ./somalierpkg/findsites
import strformat
import argparse
import algorithm
import strutils

type Site* = object
  ref_allele*: string
  alt_allele*: string
  chrom*: string
  position*: int

{.push checks: off, optimization:speed.}
proc toSite(toks: seq[string]): Site =
  result = Site()
  result.chrom = toks[0]
  result.position = parseInt(toks[1]) - 1
  result.ref_allele = toks[3]
  result.alt_allele = toks[4]

proc checkSiteRef(s:Site, fai:Fai) =
  var fa_allele = fai.get(s.chrom, s.position, s.position + s.ref_allele.len - 1).toUpperAscii
  if s.ref_allele != fa_allele:
    quit "reference base from sites file:" & s.ref_allele & " does not match that from reference: " & fa_allele
{.pop.}


proc get_sample_name(bam:Bam): string =
    var txt = newString(bam.hdr.hdr.l_text)
    copyMem(txt[0].addr, bam.hdr.hdr.text, txt.len)
    for line in txt.split("\n"):
      if line.startsWith("@RG") and "\tSM:" in line:
        var t = line.split("\tSM:")[1].split("\t")[0].strip()
        return t

    raise newException(ValueError, "[somalier] no read-group in bam file")

proc get_variant(ivcf:VCF, site:Site): Variant =
  for v in ivcf.query(&"{site.chrom}:{site.position+1}-{site.position+2}"):
    if v.start == site.position and v.REF == site.ref_allele and v.ALT[0] == site.alt_allele:
      return v.copy()

proc get_ref_alt_counts(ivcf:VCF, sites:seq[Site], fai:Fai=nil): seq[counts] =
  result = newSeq[counts](ivcf.samples.len)
  var vcf_samples = ivcf.samples
  for j, s in vcf_samples:
    result[j].sample_name = s
    result[j].sites = newSeqOfCap[allele_count](sites.len)

  var AD = newSeq[int32]()
  var n = 0

  for i, site in sites:
    var v = ivcf.get_variant(site)
    if v == nil or v.format.get("AD", AD) != Status.OK:
      AD = newSeq[int32](vcf_samples.len * 2)
    else:
      n += 1
    var mult = int(AD.len / vcf_samples.len)
    for j, s in vcf_samples:
      var ac = allele_count(nref: max(0, AD[mult * j]).uint32, nalt: max(0, AD[mult * j + 1]).uint32)

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
  # TODO: count X, Y, autosomal from sites so we can pre-allocate exactly.

  for i, site in sites:
    doAssert site.ref_allele.len == 1 and site.alt_allele.len == 1, "[somalier] can only genotype single nucleotode variants directly from alignments"
    var h = hileup(ibam, site.chrom, site.position, fai, cfg)
    checkSiteRef(site, fai)

    var ac = allele_count()

    for b in h.bases:
      if b.base.char == site.ref_allele[0]:
        ac.nref.inc
      elif b.base.char == site.alt_allele[0]:
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


proc siteOrder(a:Site, b:Site): int =
  if a.chrom == b.chrom:
    return cmp(a.position, b.position)
  return cmp(a.chrom, b.chrom)

proc readSites(path: string, fai:var Fai): seq[Site] =
  result = newSeqOfCap[Site](8192)
  var kstr = kstring_t(l:0, m:0, s:nil)
  var hf = hts_open(path.cstring, "r")

  while hts_getline(hf, cint(10), kstr.addr) > 0:
    var line  = $kstr.s
    if line[0] == '#': continue
    var sep = '\t'
    # handle ":" or tab. with ":", there is no id field.
    if sep notin line:
      sep = ':'
    var toks = line.strip().split(sep)
    if sep == ':':
      toks.insert(".", 2)

    result.add(toSite(toks))
  if len(result) > 65535:
    stderr.write_line "warning:cant use more than 65535 sites"
  sort(result, siteOrder)
  # check reference after sorting so we get much faster access.

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
    option("-d", "--out-dir", help="path to output directory")
    arg("sample_file", help="single-sample CRAM/BAM file or multi/signle-sample VCF from which to extract")

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

  var sites = readSites(opts.sites, fai)
  createDir(opts.out_dir)

  var ibam: Bam
  var ivcf: VCF

  var sample_counts: seq[counts]

  if opts.sample_file.endsWith(".vcf.gz") or opts.sample_file.endswith(".vcf.bgz") or opts.sample_file.endswith(".bcf"):
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
    "find-sites": pair(f:findsites_main, description: "create a new sites.vcf.gz file from a population VCF (this is rarely needed)."),
  }.toOrderedTable


  stderr.write_line "somalier version: " & somalierVersion
  var args = commandLineParams()

  if len(args) == 0 or not (args[0] in dispatcher):
    stderr.write_line "Commands: "
    for k, v in dispatcher:
      echo &"  {k:<13}:   {v.description}"
    if len(args) > 0 and (args[0] notin dispatcher) and args[0] notin @["-h", "-help"]:
      echo &"unknown program '{args[0]}'"
    quit ""

  dispatcher[args[0]].f()


when isMainModule:
  main()
