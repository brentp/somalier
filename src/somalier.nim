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
  ref_allele*: char
  alt_allele*: char
  chrom*: string
  position*: int

{.push checks: off, optimization:speed.}
proc toSite(toks: seq[string]): Site =
  result = Site()
  result.chrom = toks[0]
  result.position = parseInt(toks[1]) - 1
  result.ref_allele = toks[3][0]
  result.alt_allele = toks[4][0]

proc checkSiteRef(s:Site, fai:Fai) =
  var fa_allele = fai.get(s.chrom, s.position, s.position)[0].toUpperAscii
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

type counts* = object
  sample_name*: string
  sites*: seq[count]
  x_sites*: seq[count]
  y_sites*: seq[count]

proc get_ref_alt_counts(ibam:Bam, sites:seq[Site], fai:Fai): counts =

  result.sites = newSeq[relate.count](sites.len)
  result.sample_name = ibam.get_sample_name()

  var cfg = Config(MinMappingQuality: 1, ExcludeFlags:BAM_FUNMAP or BAM_FSECONDARY or BAM_FQCFAIL or BAM_FDUP)

  shallow(result.sites)
  shallow(result.x_sites)
  shallow(result.y_sites)

  var rsites = result.sites

  for i, site in sites:
    var h = hileup(ibam, site.chrom, site.position, fai, cfg)
    checkSiteRef(site, fai)

    case site.chrom:
      of ["X", "chrX"]:
        rsites = result.x_sites
      of ["Y", "chrY"]:
        rsites = result.y_sites
      else:
        rsites = result.sites

    for b in h.bases:
      if b.base.char == site.ref_allele:
        rsites[i].nref.inc
      elif b.base.char == site.alt_allele:
        rsites[i].nalt.inc
      else:
        rsites[i].nother.inc


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
    option("-t", "--threads", help="number of decompression threads to use", default="3")
    option("-d", "--out-dir", help="path to output directory")
    arg("sample_file", help="sample CRAM/BAM file from which to extract")

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
  if not open(ibam, opts.sample_file, fai=opts.fasta, index=true):
    quit "[somalier] couldn't open :" & opts.sample_file

  var cnts = ibam.get_ref_alt_counts(sites, fai)

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
    "relate": pair(f:rel_main, description: "aggregate `extract`ed information and calculate relatedness among samples"),
    "find-sites": pair(f:findsites_main, description: "create a new sites.vcf.gz file from a population VCF (this is rarely needed)"),
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
