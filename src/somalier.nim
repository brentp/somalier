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
  chrom: string
  position: int

type count* = object
  nref: int
  nalt: int

type pair = tuple[a:string, b:string, rel:float64]

proc `%`*(p:pair): JsonNode =
  return %*{"a":p.a, "b":p.b, "rel":p.rel}

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

  if ab < 0.07:
    return 0
  if ab > 0.88:
    return 2

  if ab < 0.2 or ab > 0.8: return -1 # exclude mid-range hets.

  return 1

proc count_alleles(b:Bam, site:Site): count = #{.inline.} =
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
      if over > qoff: continue
      var base = aln.base_at(qoff - over)
      if base == site.ref_allele:
        result.nref += 1
      else:
        result.nalt += 1

proc writeHelp() =
  stderr.write """
somalier [options] <bam/cram>...

version: 0.1.0

Arguments:
  <bam/cram> file(s) for samples of interest.

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


proc get_alts(bam:Bam, sites:seq[Site], nalts: ptr seq[int8], dp_stat: ptr RunningStat, min_depth:int=6): bool =
  ## count alternate alleles in a single bam at each site.
  for i, site in sites:
    var c = bam.count_alleles(site)
    dp_stat.push(c.nref + c.nalt)

    nalts[][i] = c.alts(min_depth)

proc get_bam_alts(path:string, sites:seq[Site], nalts: ptr seq[int8], dp_stat: ptr RunningStat, min_depth:int=6): bool =
  var bam: Bam
  if not open(bam, path, index=true):
    quit "couldn't open :" & $path
  result = bam.get_alts(sites, nalts, dp_stat, min_depth)
  bam.close()

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

proc checkSiteRef(s:Site, fai:var Fai) =
  var fa_allele = fai.get(s.chrom, s.position, s.position)[0]
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
  # TODO: remove s.name
  return @[s.name.split("_")[0]]

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


proc write(fh:File, sample_names: seq[string], dp_stats: seq[RunningStat], gt_counts: array[4, seq[uint16]]) =
  fh.write("#sample\tdepth_mean\tdepth_sd\tdepth_skew\tn_hom_ref\tn_het\tn_hom_alt\tn_unknown\n")
  for i, sample in sample_names:
    fh.write(&"{sample}\t{dp_stats[i].mean():.1f}\t{dp_stats[i].standard_deviation():.1f}\t{dp_stats[i].skewness()}\t{gt_counts[0][i]}\t{gt_counts[1][i]}\t{gt_counts[2][i]}\t{gt_counts[3][i]}\n")
  fh.close()

proc toj(samples: seq[string], dp_stats: seq[RunningStat], gt_counts: array[4, seq[uint16]]): string =
  result = newStringOfCap(10000)
  result.add("[")
  for i, s in samples:
    if i > 0: result.add(",\n")
    result.add($(%* {
      "sample": s,
      "depth_mean": dp_stats[i].mean(),
      "depth_std": dp_stats[i].standard_deviation(),
      "depth_skew": dp_stats[i].skewness(),
      "n_hom_ref": gt_counts[0][i],
      "n_het": gt_counts[1][i],
      "n_hom_alt": gt_counts[2][i],
      "n_unknown": gt_counts[3][i]
    }
    ))
  result.add("]")

proc main() =

  var p = initOptParser()

  var
    bv_paths = newSeq[string]()
    sites_path: string
    fai: Fai
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
      bv_paths.add(key)
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

  for i, path in bv_paths:
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
    results = newSeq[seq[int8]](n_samples)
    dp_stats = newSeq[RunningStat](n_samples)
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
    results[j] = newSeq[int8](sites.len)
    #while not preferSpawn():
    #  sleep(50)
    if responses.len > 50 and j mod 25 == 0:
      stderr.write_line "[somalier] spawning sample:", j
    responses[j] = spawn get_bam_alts(bv_paths[j], sites, results[j].addr, dp_stats[j].addr, min_depth)

  for index, fv in responses:
    blockUntil(fv)

  stderr.write_line "[somalier] collected sites from all samples"
  shallow(results)

  var t0 = cpuTime()
  var nsites = 0
  var alts = newSeq[int8](n_samples)
  # counts of hom-ref, het, hom-alt, unk
  var gt_counts : array[4, seq[uint16]]
  for i in 0..<4:
    gt_counts[i] = newSeq[uint16](n_samples)

  for rowi in 0..sites.high:
    var nun = 0
    for i in 0..<n_samples:
      alts[i] = results[i][rowi]
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

  fh_html.write(tmpl_html.replace("<INPUT_JSON>", $j).replace("<SAMPLE_JSON>", toj(sample_names, dp_stats, gt_counts)))
  fh_html.close()
  stderr.write_line("[somalier] wrote interactive HTML output to: ",  output_prefix & "html")

  fh_samples.write(sample_names, dp_stats, gt_counts)

  for rel in relatedness(final, grouped):
    fh_tsv.write_line $rel

  fh_tsv.close()
  grouped.write(output_prefix)

  stderr.write_line("[somalier] wrote groups to: ",  output_prefix & "groups.tsv")


when isMainModule:
  main()
