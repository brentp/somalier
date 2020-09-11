import strutils
import strformat
import algorithm
import hts/private/hts_concat
import hts/fai
import os

const formatVersion* = 2'u8

type Site* = object
  A_allele*: string
  B_allele*: string
  chrom*: string
  position*: int
  ## for GVCF, we don't know the actual REF/ALT allele so we can't compare to
  ## A/B_allele so instead we keep the `flip` attribute and just swap ref, atl
  ## counts if swap is true.
  flip*:bool


proc update_with_glob*(files: var seq[string]) =
  var toadd = newSeqOfCap[string](256)
  for i in 0..<min(files.len, 10):
    if files[i] == "++":
      toadd.add(files[i])
      continue
    for w in files[i].walkFiles:
      toadd.add(w)

  if files.len > 10:
    files = toadd & files[10..files.high]
  else:
    files = toadd

{.push checks: off, optimization:speed.}
proc toSite(toks: seq[string]): Site =
  result = Site()
  result.chrom = toks[0]
  result.position = parseInt(toks[1]) - 1
  if toks[3] < toks[4]:
    result.A_allele = toks[3]
    result.B_allele = toks[4]
  else:
    result.B_allele = toks[3]
    result.A_allele = toks[4]
    result.flip = true

proc checkSiteRef*(s:Site, fai:Fai) =
  var fa_allele = fai.get(s.chrom, s.position, s.position + s.A_allele.len - 1).toUpperAscii
  if s.A_allele != fa_allele and s.B_allele != fa_allele:
    quit "neither allele from sites file:" & s.A_allele & "/" & s.B_allele & " matches that from reference: " & fa_allele
{.pop.}


proc siteOrder(a:Site, b:Site): int =
  if a.chrom == b.chrom:
    return cmp(a.position, b.position)
  return cmp(a.chrom, b.chrom)

proc readSites*(path: string): seq[Site] =
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

type allele_count* = object
  nref*: uint32
  nalt*: uint32
  nother*: uint32

type counts* = object
  sample_name*: string
  sites*: seq[allele_count]
  x_sites*: seq[allele_count]
  y_sites*: seq[allele_count]

import streams

proc write_counts*(cnts: counts, sample_name: string, fname: string) =
    var s = newFileStream(fname, fmWrite)
    if s == nil:
      raise newException(IOError, &"somalier: error opening file: {fname}")
    s.write(formatVersion.uint8)
    s.write(sample_name.len.uint8)
    s.write(sample_name)
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
