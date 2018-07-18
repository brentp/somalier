{.experimental.}

import hts
import parseopt
import ospaths
import strutils
import threadpool


type Site* = object
  ref_allele: char
  chrom: string
  position: int

type count* = object
  nref: int
  nalt: int

proc ab*(c:count): float {.inline.} =
  return c.nalt.float / (c.nalt + c.nref).float

proc alts*(c:count, min_depth:int): int8 {.inline.} =
  ## give an estimate of number of alts from counts of ref and alt
  ## AB < 0.15 is called as hom-ref
  ## AB > 0.75 is hom-alt
  ## 0.15 <= AB <= 0.75 is het
  if c.nref + c.nalt < min_depth:
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
somalier rel [options] <bam/cram>...

Arguments:
  <bam/cram> file(s) for samples of interest.

Options:

  -s --sites <txt>        text file with lines of sites to use for relatedness estimation.
  -t --threads <int>        number of processors to use for parallelization.
  -f --fasta <reference>  path to reference fasta file.

  """

proc get_bam_alts(bams:seq[BAM], site:Site, nalts: var seq[int8], min_depth:int=6): bool =
  if nalts.len != 0:
    nalts.set_len(0)

  var nknown = 0
  var nref = 0

  for bam in bams:
    var c = bam.count_alleles(site)
    var calts = c.alts(min_depth)
    if calts != -1:
      nknown += 1
    if c.nref > 0:
      nref = 1
    nalts.add(calts)
  if nref == 0:
    stderr.write_line "[somalier] no reference alleles found at: ", $site


  result = nknown.float / bams.len.float >= 0.66

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

  result = nknown.float / n.float >= 0.66

proc krelated(alts: var seq[int8], ibs: var seq[uint16], n: var seq[uint16], hets: var seq[uint16], n_samples: int): int =

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

type relation_matrices = ref object
   sites_tested: int
   ibs: seq[uint16]
   n: seq[uint16]
   hets: seq[uint16]
   n_samples: int

proc relmatrix(paths:seq[string], sites:seq[Site], p: ptr relation_matrices, threads:int, idx:int): bool {.thread.} =
  result = true

  var bams = newSeq[Bam](len(paths))
  for i, p in paths:
    var b:Bam
    open(b, p, index=true, threads=threads)
    bams[i] = b

  var nalts = newSeqOfCap[int8](16)
  var rel = p[]
  if rel == nil:
    return false
  rel.sites_tested = 0

  var bad = false

  if rel.hets == nil:
    stderr.write_line "skipped:", len(sites), " ", $(rel.n == nil), " ", $(rel.ibs == nil)
    result = false
    return

  zeroMem(rel.hets[0].addr, sizeof(rel.hets[0]) * rel.hets.len)
  zeroMem(rel.n[0].addr, sizeof(rel.n[0]) * rel.n.len)
  zeroMem(rel.ibs[0].addr, sizeof(rel.ibs[0]) * rel.ibs.len)

  var n_samples = bams.len

  # TODO: allow mix of BAM/VCF
  for s in sites:
    rel.sites_tested += 1
    if not get_bam_alts(bams, s, nalts):
      continue

    discard krelated(nalts, rel.ibs, rel.n, rel.hets, n_samples)

  #p[] = rel
  #stderr.write_line "in index:", idx, " tested ", (p[]).sites_tested


iterator relatedness(r:relation_matrices, sample_names:seq[string]): relation =

  for sj in 0..<r.n_samples - 1:
    for sk in sj + 1..<r.n_samples:
      if sj == sk: quit "logic error"

      var bottom = min(r.hets[sk], r.hets[sj]).float64
      if bottom == 0:
        bottom = max(r.hets[sk], r.hets[sj]).float64
      if bottom == 0:
        # can't calculate relatedness
        bottom = -1'f64

      var grelatedness = (r.ibs[sk * r.n_samples + sj].float64 - 2 * r.ibs[sj * r.n_samples + sk].float64) / bottom

      yield relation(sample_a: sample_names[sj],
                     sample_b: sample_names[sk],
                     hets_a: r.hets[sj], hets_b: r.hets[sk],
                     ibs0: r.ibs[sj * r.n_samples + sk],
                     shared_hets: r.ibs[sk * r.n_samples + sj],
                     ibs2: r.n[sk * r.n_samples + sj],
                     n: r.n[sj * r.n_samples + sk],
                     rel: grelatedness)


proc toSite(toks: seq[string], fai:Fai): Site =
  result = Site()
  result.chrom = toks[0]
  result.position = parseInt(toks[1]) - 1
  var base = fai.get(result.chrom, result.position, result.position)
  result.ref_allele = base[0]


proc main() =

  var p = initOptParser()

  var
    bv_paths = newSeq[string]()
    sites_path: string
    fasta_path: string
    threads = 1

  for kind, key, val in p.getopt():
    case kind
    of cmdArgument:
      if key == "rel": continue
      bv_paths.add(key)
    of cmdLongOption, cmdShortOption:
      case key
      of "help", "h":
        writeHelp()
        continue
      of "threads":
        threads = parseInt(val)
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
  if threads < 1: threads = 1

  var
    fai: Fai
    sites = newSeqOfCap[Site](10000)

  if not open(fai, fasta_path):
    quit "couldn't open fasta with fai:" & fasta_path

  for line in sites_path.lines:
    var toks = line.strip().split(":")
    #if toks[0] != last_chrom:
    #  #echo "getting chrom:", toks[0]
    #  last_chrom = toks[0]
    #  cseq = fai.get(last_chrom)

    sites.add(toSite(toks, fai))

  if len(sites) > 65535:
    quit "cant use more than 65535 sites"

  var sample_names = newSeq[string](len(bv_paths))

  for i, path in bv_paths:
    #var bam:Bam
    #open(bam, path, index=true, threads=2)
    #bams[i] = bam
    var s = splitFile(path)
    sample_names[i] = s.name
  var
    n_samples = sample_names.len
    batch_size = 8 #(len(sites) / threads).int


  if threads > len(sites):
    quit "can't use more threads than sites"
  if threads < 2:
    batch_size = len(sites)
  if threads * batch_size > sites.len:
    stderr.write_line "[somalier] setting to lower number of threads to avoid tiny batches of work"
    threads = (sites.len / batch_size).int
  var jobs = (sites.len / batch_size).int

  var
    results = newSeq[relation_matrices](threads)
    responses = newSeq[FlowVarBase](threads)


  stderr.write_line "[somalier] batch-size:", batch_size, " sites:", len(sites), " threads:" & $threads, " jobs:", $jobs

  # aggregated from all threads.
  var final = relation_matrices(ibs: newSeq[uint16](n_samples * n_samples),
                              n: newSeq[uint16](n_samples * n_samples),
                              hets: newSeq[uint16](n_samples),
                              n_samples:n_samples)

  for j in 0..<responses.len:
    results[j] = relation_matrices()
    results[j].ibs = newSeq[uint16](n_samples * n_samples)
    results[j].n = newSeq[uint16](n_samples * n_samples)
    results[j].hets = newSeq[uint16](n_samples)
    results[j].n_samples = n_samples
    results[j].sites_tested = 0
    responses[j] = spawn relmatrix(bv_paths, sites[(j * batch_size)..<((j + 1) * batch_size)], results[j].addr, 1, j)


  var jobi = responses.len
  while results.len != 0:
    var index = awaitAny(responses)
    if index == -1:
      quit "[somalier] got unexpected value from await"
    var relm = results[index]

    if final.n_samples != relm.n_samples:
      quit "[somalier] got differing numbers of samples"

    #stderr.write_line "got ", $relm.sites_tested, " for index ", $index
    final.sites_tested += relm.sites_tested

    for i, v in relm.n:
      final.n[i] += v
      if relm.sites_tested == 0 and v > 0'u16:
        quit "WTF n"
    for i, v in relm.hets:
      final.hets[i] += v
      if relm.sites_tested == 0 and v > 0'u16:
        quit "WTF hets"
    for i, v in relm.ibs:
      final.ibs[i] += v
      if relm.sites_tested == 0 and v > 0'u16:
        quit "WTF ibs"

    # send next job
    if (jobi * batch_size) < sites.len:
      var
        imin = jobi * batch_size
        imax = min(sites.len, (jobi + 1) * batch_size)
      #echo "sending ", imin, "..", imax, " of ", sites.len, " to index:", index
      responses[index] = spawn relmatrix(bv_paths, sites[imin..<imax], results[index].addr, 2, jobi)
    else:
      results.del(index)
      responses.del(index)

    jobi += 1
    
  stderr.write_line "[somalier] sites tested:", $final.sites_tested
  for rel in relatedness(final, sample_names):
    echo rel


when isMainModule:
  GC_disableMarkAndSweep()
  main()
