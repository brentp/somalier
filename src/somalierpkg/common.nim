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

type SiteAF* = tuple[site: Site, af: float32]

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

proc siteAFOrder(a: SiteAF, b: SiteAF): int =
  siteOrder(a.site, b.site)

proc isAutosomal*(site: Site): bool =
  site.chrom notin ["X", "chrX", "NC_000023.10", "NC_000023.11",
                    "Y", "chrY", "NC_000024.9", "NC_000024.10"]

proc stored_alt_af*(site: Site, vcf_af: float32): float32 =
  ## Somalier stores counts in alphabetical A/B allele order, so sites whose
  ## REF/ALT order was flipped during extraction need the population AF flipped
  ## here as well to stay aligned with the stored B allele counts.
  if site.flip: 1'f32 - vcf_af else: vcf_af

proc parse_info_af(info: string, path: string): float32 =
  for field in info.split(';'):
    if not field.startsWith("AF="):
      continue
    let value = field[3 .. field.high]
    if value.len == 0:
      break
    try:
      return parseFloat(value.split(',')[0]).float32
    except ValueError:
      quit "malformed AF in sites file: " & path

  quit "missing AF in sites file: " & path

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

proc readSitesWithAF*(path: string): tuple[sites: seq[Site], pop_afs: seq[float32]] =
  var pairs = newSeqOfCap[SiteAF](8192)
  var kstr = kstring_t(l:0, m:0, s:nil)
  var hf = hts_open(path.cstring, "r")

  while hts_getline(hf, cint(10), kstr.addr) > 0:
    var line = $kstr.s
    if line.len == 0 or line[0] == '#':
      continue
    if '\t' notin line:
      quit "sites with AF must be tab-delimited VCF/BCF-like records: " & path

    let toks = line.strip().split('\t')
    if toks.len < 8:
      quit "expected at least 8 columns in sites file: " & path

    let site = toSite(toks)
    if not site.isAutosomal:
      continue
    # Parse AF before sorting, then sort the (site, af) pairs together so the
    # returned AF vector stays in exactly the same order as the autosomal sketch.
    pairs.add((site: site, af: stored_alt_af(site, parse_info_af(toks[7], path))))

  if pairs.len > 65535:
    stderr.write_line "warning:cant use more than 65535 sites"
  sort(pairs, siteAFOrder)
  result.sites = newSeq[Site](pairs.len)
  result.pop_afs = newSeq[float32](pairs.len)
  for i, pair in pairs:
    result.sites[i] = pair.site
    result.pop_afs[i] = pair.af

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
