import parseopt
import ospaths
import hts
import strutils


type Site* = object
  ref_allele: char
  chrom: string
  position: int

type count* = object
  nref: int
  nalt: int

proc ab*(c:count): float {.inline.} =
  return c.nalt.float / (c.nalt + c.nref).float

proc alts*(c:count): int8 {.inline.} =
  ## give an estimate of number of alts from counts of ref and alt
  ## AB < 0.15 is called as hom-ref
  ## AB > 0.75 is hom-alt
  ## 0.15 <= AB <= 0.75 is het
  if c.nref + c.nalt < 5:
    return -1
  if c.nalt == 0:
    return 0

  var ab = c.ab

  if ab < 0.05:
    return 0
  if ab > 0.9:
    return 2

  if ab < 0.25 or ab > 0.75: return -1 # exclude mid-range hets.

  return 1

proc count_alleles(b:Bam, site:Site): count =
  for aln in b.query(site.chrom, site.position, site.position + 1):
    var off = aln.start
    var qoff = 0
    var roff_only = 0
    for event in aln.cigar:
      var cons = event.consumes
      if cons.query:
        qoff += event.len
      if cons.reference:
        off += event.len
        if not cons.query:
          roff_only += event.len
      if off <= site.position: continue

      var over = off - site.position - roff_only
      var base = aln.base_at(qoff - over)
      if base == site.ref_allele:
        result.nref += 1
      else:
        result.nalt += 1

proc writeHelp() =
  stderr.write """
somalier vcf [options] <vcf/bcf>...

Arguments:
  <vcf/bcf> file(s) containing variant calls

Options:

  -s --sites <txt>    optional text file with lines of sites to use for relatedness estimation.

  """

proc get_bam_alts(bams:seq[BAM], site:Site, nalts: var seq[int8]): bool =
  if nalts.len != 0:
    nalts.set_len(0)

  var nknown = 0
  var nref = 0

  for bam in bams:
    var c = bam.count_alleles(site)
    if c.alts != -1:
      nknown += 1
    if c.nref > 0:
      nref = 1
    nalts.add(c.alts)
  if nref == 0:
    stderr.write_line "[somalier] no reference alleles found at: ", $site


  result = nknown.float / bams.len.float >= 0.8

proc get_alts(vcfs:seq[VCF], site:string, nalts: var seq[int8], cache: var seq[int32]): bool =
  if nalts.len != 0:
    nalts.set_len(0)

  var sp = site.split(":")
  var p1 = parseInt(sp[1])
  var q = sp[0] & ":" & $p1 & "-" & $(p1 + 1)

  var nknown: int
  var n:int

  for vcf in vcfs:
    var found = 0
    n += vcf.n_samples
    for v in vcf.query(q):
      if v.REF != sp[2]: continue
      var alts = v.format.genotypes(cache).alts
      for alt in alts:
        if alt != -1:
          found += 1
      nalts.add(alts)

    nknown += found

  result = nknown.float / n.float >= 0.8

proc krelated(alts: seq[int8], ibs: var seq[uint16], n: var seq[uint16], hets: var seq[uint16], n_samples: int): int =

  if alts[n_samples - 1] == 1:
    hets[n_samples-1] += 1

  var is_het: bool
  var aj, ak: int8
  var nused = 0

  for j in 0..<(n_samples-1):
    aj = alts[j]
    if aj == -1: continue
    is_het = (aj == 1)

    if is_het:
      hets[j] += 1
    nused += 1

    for k in j+1..<n_samples:
      ak = alts[k]
      if ak == -1: continue
      n[j * n_samples + k] += 1
      if is_het:
        # shared hets
        if ak == 1:
          ibs[k * n_samples + j] += 1
      else:
        # ibs0
        if aj != ak and aj + ak == 2:
          ibs[j * n_samples + k] += 1
      # ibs2
      if aj == ak:
        n[k * n_samples + j] += 1
  return nused

type relation = object
  sample_a: string
  sample_b: string
  hets_a: uint16
  hets_b: uint16
  shared_hets: uint16
  rel: float64
  ibs0: uint16
  ibs2: uint16
  n: uint16

iterator site_relatedness(bams:seq[Bam], sample_names: seq[string], sites:seq[Site]): relation =

  var nalts = newSeqOfCap[int8](16)

  var n_samples = bams.len

  var ibs = newSeq[uint16](n_samples * n_samples)
  var n = newSeq[uint16](n_samples * n_samples)
  var hets = newSeq[uint16](n_samples)

  for i, s in sites:
    echo i, "/", len(sites), " ", s
    if not get_bam_alts(bams, s, nalts):
      continue

    discard krelated(nalts, ibs, n, hets, n_samples)

  for sj in 0..<n_samples - 1:
    for sk in sj + 1..<n_samples:
      if sj == sk: quit "logic error"

      var bottom = min(hets[sk], hets[sj]).float64
      if bottom == 0:
        bottom = max(hets[sk], hets[sj]).float64
      if bottom == 0:
        # can't calculate relatedness
        bottom = -1'f64

      var relatedness = (ibs[sk * n_samples + sj].float64 - 2 * ibs[sj * n_samples + sk].float64) / bottom

      yield relation(sample_a: sample_names[sj],
                     sample_b: sample_names[sk],
                     hets_a: hets[sj], hets_b: hets[sk],
                     ibs0: ibs[sj * n_samples + sk],
                     shared_hets: ibs[sk * n_samples + sj],
                     ibs2: n[sk * n_samples + sj],
                     n: n[sj * n_samples + sk],
                     rel: relatedness)



proc toSite(toks: seq[string], rseq:var string): Site =
  result = Site()
  result.chrom = toks[0]
  result.position = parseInt(toks[1]) - 1
  result.ref_allele = rseq[result.position]


proc main() =

  var p = initOptParser()

  var
    vcf_paths = newSeq[string]()
    sites_path: string
    fasta_path: string

  for kind, key, val in p.getopt():
    case kind
    of cmdArgument:
      vcf_paths.add(key)
    of cmdLongOption, cmdShortOption:
      case key
      of "help", "h":
        writeHelp()
      of "vcf":
        continue
      of "sites", "s":
        sites_path = val
      of "fasta", "f":
        fasta_path = val
      else:
        writeHelp()
    of cmdEnd:
      assert(false)
  if sites_path == nil or sites_path == "":
    echo "must set sites path"
    writeHelp()
    quit(2)
  if fasta_path == nil or fasta_path == "":
    echo "must set fasta path"
    writeHelp()
    quit(2)

  

  var
    fai: Fai
    cseq: string
    sites = newSeqOfCap[Site](10000)
    last_chrom: string = ""

  if not open(fai, fasta_path):
    quit "couldn't open fasta with fai:" & fasta_path

  for line in sites_path.lines:
    var toks = line.strip().split(":")
    if toks[0] != last_chrom:
      echo "getting chrom:", toks[0]
      last_chrom = toks[0]
      cseq = fai.get(last_chrom)

    sites.add(toSite(toks, cseq))

  if len(sites) > 65535:
    quit "cant use more than 65535 sites"

  var vcfs = newSeq[VCF](len(vcf_paths))
  var bams = newSeq[Bam](len(vcf_paths))
  var sample_names = newSeq[string](len(vcf_paths))

  for i, path in vcf_paths:
    var bam:Bam
    open(bam, path, index=true, threads=2)
    bams[i] = bam
    var s = splitFile(path)
    sample_names[i] = s.name

  for rel in site_relatedness(bams, sample_names, sites):
    echo rel

when isMainModule:
  main()
