import strutils
import algorithm
import hts/private/hts_concat
import hts/fai

type Site* = object
  A_allele*: string
  B_allele*: string
  chrom*: string
  position*: int

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
