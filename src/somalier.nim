{.experimental.}

import os
import hts
import sets
import math
import json
import stats
import times
import bpbiopkg/pedfile
import strformat
import ./parseopt3
import algorithm
import strutils
import threadpool
import ./results_html

type Site* = object
  ref_allele: char
  alt_allele: char
  chrom: string
  position: int

type count* = object
  nref: uint32
  nalt: uint32
  nother: uint32

type pair = tuple[a:string, b:string, rel:float64]

proc `%`*(p:pair): JsonNode =
  return %*{"a":p.a, "b":p.b, "rel":p.rel}

template proportion_other(c:count): float =
  if c.nother == 0: 0'f else: c.nother.float / (c.nother + c.nref + c.nalt).float

proc ab*(c:count, min_depth:int): float {.inline.} =
  if c.proportion_other > 0.04: return -1
  if int(c.nref + c.nalt) < min_depth:
    return -1
  if c.nalt == 0:
    return 0
  result = c.nalt.float / (c.nalt + c.nref).float

proc alts(ab:float): int8 {.inline.} =
  if ab < 0: return -1
  if ab < 0.07:
    return 0
  if ab > 0.88:
    return 2
  if ab < 0.2 or ab > 0.8: return -1 # exclude mid-range hets.
  return 1

proc count_alleles(b:Bam, site:Site): count {.inline.} =
  for aln in b.query(site.chrom, site.position, site.position + 1):
    if aln.mapping_quality < 10: continue
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

      if not (cons.query and cons.reference): continue

      var over = off - site.position - roff_only
      if over > qoff: break
      if over < 0: continue
      doAssert qoff - over >= 0
      var base = aln.base_at(qoff - over)
      if base == site.ref_allele:
        result.nref += 1
      elif base == site.alt_allele:
        result.nalt += 1
      else:
        #stderr.write_line $event, " -> over:", over, " -> ", site, " -> ", aln.tostring
        result.nother += 1
    doAssert aln.start >= 0

proc writeHelp() =
  stderr.write """
somalier [options] <bam/cram/list>...

version: 0.1.4

Arguments:
  <bam/cram/list> file(s) for samples of interest. or a file ending in ".list" where each line
  contains the path to the (possibly remote) alignment file in column 1 and optionally the index
  in column 2. (tab, space, or comma delimited).

Options:

  -s --sites <vcf>        vcf file with lines of sites to use for relatedness estimation.
  -t --threads <int>      optional number of processors to use for parallelization.
  -d --min-depth <int>    only consider sites with at least this depth [default: 7].
  -f --fasta <reference>  path to reference fasta file.
  -g --groups <path>      optional path to expected groups of samples (e.g. tumor normal pairs).
                          specified as comma-separated groups per line e.g.:
                            normal1,tumor1a,tumor1b
                            normal2,tumor2a
  -p --ped <path>         optional path to a ped/fam file indicating the expected relationships
                          among samples.
  -o --output <prefix>    output prefix for results.

  """

type Stat4 = object
  dp: RunningStat
  gtdp: RunningStat # depth of genotyped sites
  un: RunningStat
  ab: RunningStat

type pathWithIndex = object
  ## this allows us to support paths where the index and the bam are in different locations
  path: string
  indexPath: string

const arraySize = 1024

type counter = object
  counts: array[arraySize, int]

proc push(c:var counter, x:float) {.inline.} =
  doAssert 0 <= x and x <= 0.5
  var idx = min(c.counts.high, int(x * (arraySize * 2.0) + 0.5))
  c.counts[idx].inc

proc median(c:counter): float {.inline.} =
  var n = sum(c.counts)
  if n < 1: return 0
  var mid = int(n.float / 2.0 + 0.5)
  var fsum = 0
  for i, v in c.counts:
    fsum += v
    if fsum >= mid:
      return i.float / (arraySize * 2.0)

proc get_abs(ibam:Bam, sites:seq[Site], ab: ptr seq[float32], stat: ptr Stat4, min_depth:int=6): bool =
  ## count alternate alleles in a single bam at each site.
  for i, site in sites:
    var c = ibam.count_alleles(site)
    stat.dp.push(int(c.nref + c.nalt))
    if c.nref > 0'u32 or c.nalt > 0'u32 or c.nother > 0'u32:
      stat.un.push(c.nother.float64 / float64(c.nref + c.nalt + c.nother))
    if c.nref.float > min_depth / 2 and c.nalt.float > min_depth / 2:
      stat.ab.push(c.ab(min_depth))

    ab[][i] = c.ab(min_depth)
    if ab[][i] != -1:
      stat.gtdp.push(int(c.nref + c.nalt))

const error_rate = 2e-3

proc estimate_contamination(self_abs: seq[float32], other_abs: seq[float32]): (float, int) =
  ## estimate contamination of self, by other.
  var sites_used = 0
  var c: counter
  for i, a in self_abs:
    if a == -1: continue
    var b = other_abs[i]
    if b == -1: continue
    if abs(a - b) < 0.015: continue

    sites_used += 1
    # aa: 0.5, bb: 0 -> 2
    # aa: 1,   bb: 0 -> 1
    # aa: 0.01, bb: 0.03 -> 2
    if a > 0.25 and a < 0.75: continue

    var scaler = min(2'f32, 1'f32 / abs(a - b).float32)
    # a: 0.01, b: 0.0001 then b can't contribute to a
    if a < 0.5 and b < a:
      scaler = 0

    # a: 0.99, b: 0.999 then reads can't come from b
    elif a > 0.5 and b > a:
      scaler = 0

    var ax = if a > 0.5: 1-a else: a
    if ax > 0.25:
      # remove cases where direction of support does not match
      # e.g.:
      # self:0.4744310677051544 other:0.9432255029678345 ax:0.02556893229484558 scaler:2.0 evidence for contamination of:0.05113786458969116
      if a < 0.5 and b > 0.5:
        scaler = 0.0 # no way that change in AF can come from b
      ax = 0.5 - ax

    var contam = scaler * ax

    #echo "self:", a, " other:", b, " ax:", ax, " scaler:", scaler, " evidence for contamination of:", contam

    c.push(contam) #scaler * (if a < 0.5: a else: 1 - a))
  #echo sum(c.counts), " ", c.counts[0..<min(arraySize, 100)]
  return (c.median, sum(c.counts))

proc get_bam_abs(pwi:pathWithIndex, fai:string, sites:seq[Site], ab: ptr seq[float32], stat: ptr Stat4, min_depth:int): bool {.thread.} =
  var ibam: Bam
  if not open(ibam, pwi.path, fai=fai):
    quit "couldn't open :" & $pwi.path
  ibam.load_index(pwi.indexPath)
  discard ibam.set_option(FormatOption.CRAM_OPT_REQUIRED_FIELDS, 8191 - SAM_QUAL.int - SAM_QNAME.int)
  result = ibam.get_abs(sites, ab, stat, min_depth)
  doAssert ibam.idx != nil
  ibam.close()

proc get_depths(v:Variant, cache: var seq[int32]): seq[int32] =
  for c in cache.mitems:
    c = 0
  if v.format.get("AD", cache) == Status.OK:
    result = newSeq[int32](v.n_samples)
    for i in 0..<v.n_samples:
      result[i] = cache[2*i] + cache[2*i+1]
    return

  if v.format.get("DP", cache) == Status.OK:
    result = newSeq[int32](v.n_samples)
    copyMem(result[0].addr, cache[0].addr, sizeof(cache[0]) * cache.len)
    return


proc get_vcf_alts(vcfs:seq[VCF], site:Site, nalts: var seq[int8], cache: var seq[int32], min_depth:int): int =
  if len(vcfs) == 0:
    return 0

  var n:int

  for vcf in vcfs:
    var found = 0
    n += vcf.n_samples
    {.gcsafe.}:
      for v in vcf.query(site.chrom & ":" & $(site.position + 1) & "-" & $(site.position + 2)):
        if v.REF[0] != site.ref_allele: continue
        var alts = v.format.genotypes(cache).alts
        var dps = get_depths(v, cache)
        for k, alt in alts.mpairs:
          if dps.len != 0 and dps[k] < min_depth:
            alt = -1
          if alt != -1:
            found += 1
        nalts.add(alts)

    result += found


proc krelated(alts: var seq[int8], ibs: var seq[uint16], n: var seq[uint16], hets: var seq[uint16], homs: var seq[uint16], shared_hom_alts: var seq[uint16], n_samples: int): int {.inline.} =

  if alts[n_samples - 1] == 1:
    hets[n_samples-1] += 1
  elif alts[n_samples - 1] == 2:
    homs[n_samples-1] += 1

  var is_het: bool
  var aj, ak: int8
  var nused = 0

  for j in 0..<(n_samples-1):
    aj = alts[j]
    if aj == -1: continue
    is_het = (aj == 1)

    if is_het:
      hets[j] += 1
    elif aj == 2:
      homs[j] += 1

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
      if aj == ak: #and not is_het:
        n[k * n_samples + j] += 1
        if aj == 2:
          shared_hom_alts[j * n_samples + k] += 1
  return nused

type relation = object
  sample_a: string
  sample_b: string
  hets_a: uint16
  hets_b: uint16
  hom_alts_a: uint16
  hom_alts_b: uint16
  shared_hom_alts: uint16
  shared_hets: uint16
  ibs0: uint16
  ibs2: uint16
  n: uint16

proc hom_alt_concordance(r: relation): float64 {.inline.} =
  return (r.shared_hom_alts.float64 - 2 * r.ibs0.float64) / min(r.hom_alts_a, r.hom_alts_b).float64

proc rel(r:relation): float64 {.inline.} =
  return (r.shared_hets.float64 - 2 * r.ibs0.float64) / min(r.hets_a, r.hets_b).float64

const header = "$sample_a\t$sample_b\t$relatedness\t$hom_concordance\t$hets_a\t$hets_b\t$shared_hets\t$hom_alts_a\t$hom_alts_b\t$shared_hom_alts\t$ibs0\t$ibs2\t$n"
proc `$`(r:relation): string =
  return header % [
         "sample_a", r.sample_a, "sample_b", r.sample_b,
         "relatedness", formatFloat(r.rel, ffDecimal, precision=3),
         "hom_concordance", formatFloat(r.hom_alt_concordance, ffDecimal, precision=3),
         "hets_a", $r.hets_a, "hets_b", $r.hets_b,
         "shared_hets", $r.shared_hets, "hom_alts_a", $r.hom_alts_a, "hom_alts_b", $r.hom_alts_b, "shared_hom_alts", $r.shared_hom_alts, "ibs0", $r.ibs0, "ibs2", $r.ibs2, "n", $r.n]

type relation_matrices = object
   sites_tested: int
   ibs: seq[uint16]
   n: seq[uint16]
   hets: seq[uint16]
   homs: seq[uint16]
   shared_hom_alts: seq[uint16]
   samples: seq[string]

proc n_samples(r: relation_matrices): int {.inline.} =
  return r.samples.len

proc bam_like(path:string): bool {.inline.} =
    return path.endsWith(".bam") or path.endsWith(".cram")


iterator relatedness(r:relation_matrices, grouped: var seq[pair]): relation =
  var sample_names = r.samples

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
      if grelatedness > 0.2:
        grouped.add((sample_names[sj], sample_names[sk], grelatedness))
      yield relation(sample_a: sample_names[sj],
                     sample_b: sample_names[sk],
                     hets_a: r.hets[sj], hets_b: r.hets[sk],
                     hom_alts_a: r.homs[sj], hom_alts_b: r.homs[sk],
                     ibs0: r.ibs[sj * r.n_samples + sk],
                     shared_hets: r.ibs[sk * r.n_samples + sj],
                     shared_hom_alts: r.shared_hom_alts[sj * r.n_samples + sk],
                     ibs2: r.n[sk * r.n_samples + sj],
                     n: r.n[sj * r.n_samples + sk])


{.push checks: off, optimization:speed.}
proc toSite(toks: seq[string]): Site =
  result = Site()
  result.chrom = toks[0]
  result.position = parseInt(toks[1]) - 1
  result.ref_allele = toks[3][0]
  result.alt_allele = toks[4][0]

proc checkSiteRef(s:Site, fai:var Fai) =
  var fa_allele = fai.get(s.chrom, s.position, s.position)[0].toUpperAscii
  if s.ref_allele != fa_allele:
    quit "reference base from sites file:" & s.ref_allele & " does not match that from reference: " & fa_allele
{.pop.}

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
    if line.count(sep) == 0:
      sep = ':'
    var toks = line.strip().split(sep)
    if sep == ':':
      toks.insert(".", 2)

    result.add(toSite(toks))
  if len(result) > 65535:
    stderr.write_line "warning:cant use more than 65535 sites"
  sort(result, siteOrder)
  # check reference after sorting so we get much faster access.
  for i, r in result:
    if i mod 10000 == 0 and i > 0:
      stderr.write_line "[somalier] checked reference for " & $i & " sites"

    checkSiteRef(r, fai)
  fai = nil


proc `%`*(v:uint16): JsonNode =
  new(result)
  result.kind = JInt
  result.num = v.int64

proc readGroups(path:string): seq[pair] =
  result = newSeq[pair]()
  if path == "":
    return

  # expand out a,b,c to a,b, a,c, b,c
  for line in path.lines:
    var row = line.strip().split(",")
    var rel = 1.0
    if '\t' in row[row.high]:
      var tmp = row[row.high].split('\t')
      doAssert tmp.len == 2
      row[row.high] = tmp[0]
      rel = parseFloat(tmp[1])
    for i, x in row[0..<row.high]:
      for j, y in row[(i+1)..row.high]:
        if x < y:
          result.add((x, y, rel))
        else:
          result.add((y, x, rel))

proc get_sample_names(path: string): seq[string] =
  if path.bam_like:
    var bam: Bam
    open(bam, path)
    var txt = newString(bam.hdr.hdr.l_text)
    copyMem(txt[0].addr, bam.hdr.hdr.text, txt.len)
    for line in txt.split("\n"):
      if line.startsWith("@RG") and "\tSM:" in line:
        var t = line.split("\tSM:")[1].split("\t")[0].strip()
        return @[t]

  elif path.endsWith("vcf.gz") or path.endswith(".bcf") or path.endsWith(".bcf.gz") or path.endsWith("vcf.bgz"):
    var vcf: VCF
    if not open(vcf, path):
      quit "could not open " & $path
    return vcf.samples

  stderr.write_line "[somalier] warning couldn't find samples for " & path & " using file names to guess."
  var s = splitFile(path)
  return @[s.name]

proc write(grouped: seq[pair], output_prefix:string) =
  if len(grouped) == 0: return
  var fh_groups:File
  if not open(fh_groups,  output_prefix & "groups.tsv", fmWrite):
    quit "couldn't open groups file."
  for grp in grouped:
    fh_groups.write(&"{grp.a},{grp.b}\t{grp.rel}\n")
  fh_groups.close()

proc add_ped_samples(grouped: var seq[pair], samples:seq[Sample], sample_names:seq[string]) =
  ## samples were parsed from ped. we iterate over them and add any pair where both samples are in sample_names
  if samples.len == 0: return
  var ss = initSet[string]()
  for s in sample_names: ss.incl(s) # use a set for better lookup.
  for i, sampleA in samples[0..<samples.high]:
    if sampleA.id notin ss: continue
    for j, sampleB in samples[i + 1..samples.high]:
       if sampleB.id notin ss: continue
       var rel = relatedness(sampleA, sampleB, samples)
       if rel <= 0: continue
       if sampleA.id < sampleB.id:
         grouped.add((sampleA.id, sampleB.id, rel))
       else:
         grouped.add((sampleB.id, sampleA.id, rel))


proc write(fh:File, sample_names: seq[string], stats: seq[Stat4], gt_counts: array[4, seq[uint16]]) =
  fh.write("#sample\tgt_depth_mean\tgt_depth_sd\tgt_depth_skew\tdepth_mean\tdepth_sd\tdepth_skew\tab_mean\tab_std\tn_hom_ref\tn_het\tn_hom_alt\tn_unknown\n")
  for i, sample in sample_names:
    fh.write(&"{sample}\t")
    fh.write(&"{stats[i].gtdp.mean():.1f}\t{stats[i].gtdp.standard_deviation():.1f}\t{stats[i].gtdp.skewness()}\t")
    fh.write(&"{stats[i].dp.mean():.1f}\t{stats[i].dp.standard_deviation():.1f}\t{stats[i].dp.skewness()}\t")
    fh.write(&"{stats[i].ab.mean():.1f}\t{stats[i].ab.standard_deviation():.1f}\t{gt_counts[0][i]}\t{gt_counts[1][i]}\t{gt_counts[2][i]}\t{gt_counts[3][i]}\n")
  fh.close()

proc toj(samples: seq[string], stats: seq[Stat4], gt_counts: array[4, seq[uint16]]): string =
  result = newStringOfCap(10000)
  result.add("[")
  for i, s in samples:
    if i > 0: result.add(",\n")
    result.add($(%* {
      "sample": s,

      "gt_depth_mean": stats[i].gtdp.mean(),
      "gt_depth_std": stats[i].gtdp.standard_deviation(),
      "gt_depth_skew": stats[i].gtdp.skewness(),

      "depth_mean": stats[i].dp.mean(),
      "depth_std": stats[i].dp.standard_deviation(),
      "depth_skew": stats[i].dp.skewness(),

      "ab_mean": stats[i].ab.mean(),
      "ab_std": stats[i].ab.standard_deviation(),
      "ab_skew": stats[i].ab.skewness(),
      "pct_other_alleles": 100.0 * stats[i].un.mean,
      "n_hom_ref": gt_counts[0][i],
      "n_het": gt_counts[1][i],
      "n_hom_alt": gt_counts[2][i],
      "n_unknown": gt_counts[3][i],
      "n_known": gt_counts[0][i] + gt_counts[1][i] + gt_counts[2][i]
    }
    ))
  result.add("]")

proc add_pl(bv: var seq[pathWithIndex], arg: string) =
  ## normalize arguments between *.bam and a .list with sample, index pairs
  if arg.endswith(".list"):
    for l in arg.lines:
      var toks = l.split(seps={',', '\t', ' '})
      doAssert toks.len in {1, 2}
      if toks.len == 1: toks.add("")
      bv.add(pathWithIndex(path: toks[0], indexPath: toks[1]))
  else:
    bv.add(pathWithIndex(path:arg))

proc main() =

  var p = initOptParser()

  var
    bv_paths = newSeq[pathWithIndex]()
    sites_path: string
    fai: Fai
    fasta_path: string
    samples: seq[Sample]
    min_depth = 7
    groups: seq[pair]
    orig_groups: seq[pair]
    output_prefix: string = "somalier."
    threads = 1

  for kind, key, val in p.getopt():
    case kind
    of cmdArgument:
      if key == "rel": continue
      bv_paths.add_pl(key)
    of cmdLongOption, cmdShortOption:
      case key
      of "help", "h":
        writeHelp()
        quit(0)
      of "d", "min-depth":
        min_depth = parseInt(val)
      of "threads", "t":
        threads = parseInt(val)
      of "ped", "p":
        samples = parse_ped(val)
      of "output", "o":
        output_prefix = val.strip(chars={'.'}) & "." & "somalier."
      of "sites", "s":
        sites_path = val
      of "fasta", "f":
        if not open(fai, val):
          quit "couldn't open fasta with fai:" & val
        fasta_path = val
      of "groups", "g":
        orig_groups = readGroups(val)
      else:
        writeHelp()
    of cmdEnd:
      assert(false)
  if sites_path == "":
    echo "must set sites path"
    writeHelp()
    quit(2)
  if fai == nil:
    echo "must set fasta path"
    writeHelp()
    quit(2)
  if threads < 1: threads = 1

  var sites = readSites(sites_path, fai)

  ## need to track samples names from bams first, then vcfs since
  ## thats the order for the alts array.
  var sample_names = newSeqOfCap[string](len(bv_paths))
  var non_bam_sample_names = newSeqOfCap[string](len(bv_paths))

  for i, pw in bv_paths:
    var path = pw.path
    if path.bam_like:
      sample_names.add(get_sample_names(path))
    else:
      non_bam_sample_names.add(get_sample_names(path))

  sample_names.add(non_bam_sample_names)
  groups.add_ped_samples(samples, sample_names)
  # add orig_groups last so the -g takes precedence over -p
  groups.add(orig_groups)

  var
    n_samples = sample_names.len
  if threads < 2:
    threads = 1

  if threads > n_samples:
    threads = n_samples

  var
    ab_results = newSeq[seq[float32]](n_samples)
    stats = newSeq[Stat4](n_samples)
    responses = newSeq[FlowVarBase](n_samples)

  stderr.write_line "[somalier] sites:", len(sites), " threads:" & $threads
  #setMinPoolSize(threads)
  setMaxPoolSize(threads)

  # aggregated from all threads.
  var final = relation_matrices(ibs: newSeq[uint16](n_samples * n_samples),
                              n: newSeq[uint16](n_samples * n_samples),
                              shared_hom_alts: newSeq[uint16](n_samples * n_samples),
                              hets: newSeq[uint16](n_samples),
                              homs: newSeq[uint16](n_samples),
                              samples: sample_names)

  for j in 0..<responses.len:
    ab_results[j] = newSeq[float32](sites.len)
    if responses.len > 50 and j mod 25 == 0:
      stderr.write_line "[somalier] spawning sample:", j
    responses[j] = spawn get_bam_abs(bv_paths[j], fasta_path, sites, ab_results[j].addr, stats[j].addr, min_depth)

  for index, fv in responses:
    blockUntil(fv)

  stderr.write_line "[somalier] collected sites from all samples"
  shallow(ab_results)

  var t0 = cpuTime()
  var nsites = 0
  var alts = newSeq[int8](n_samples)
  # counts of hom-ref, het, hom-alt, unk
  var gt_counts : array[4, seq[uint16]]
  for i in 0..<4:
    gt_counts[i] = newSeq[uint16](n_samples)

  for i in 1..<ab_results.len:
    for o in 0..<i:
      var a = ab_results[i]
      var b = ab_results[o]
      echo sample_names[o], " vs ", sample_names[i], " =>", estimate_contamination(a, b)

  for rowi in 0..sites.high:
    var nun = 0
    for i in 0..<n_samples:
      alts[i] = ab_results[i][rowi].alts
      if alts[i] == -1:
        nun.inc
        gt_counts[3][i].inc
      else:
        gt_counts[alts[i].int][i].inc

    if nun.float64 / n_samples.float64 > 0.6: continue
    nsites += 1

    discard krelated(alts, final.ibs, final.n, final.hets, final.homs, final.shared_hom_alts, n_samples)

  stderr.write_line &"time to calculate relatedness on {nsites} usable sites: {cpuTime() - t0:.3f}"
  var
    fh_tsv:File
    fh_samples:File
    fh_html:File
    grouped: seq[pair]

  if not open(fh_tsv, output_prefix & "pairs.tsv", fmWrite):
    quit "couldn't open output file"
  if not open(fh_samples, output_prefix & "samples.tsv", fmWrite):
    quit "couldn't open output file"
  if not open(fh_html, output_prefix & "html", fmWrite):
    quit "couldn't open html output file"

  fh_tsv.write_line '#', header.replace("$", "")

  var j = % final
  j["expected-relatedness"] = %* groups

  fh_html.write(tmpl_html.replace("<INPUT_JSON>", $j).replace("<SAMPLE_JSON>", toj(sample_names, stats, gt_counts)))
  fh_html.close()
  stderr.write_line("[somalier] wrote interactive HTML output to: ",  output_prefix & "html")

  fh_samples.write(sample_names, stats, gt_counts)

  for rel in relatedness(final, grouped):
    fh_tsv.write_line $rel

  fh_tsv.close()
  grouped.write(output_prefix)

  stderr.write_line("[somalier] wrote groups to: ",  output_prefix & "groups.tsv")


when isMainModule:
  main()
