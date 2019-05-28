import os
import strformat
import times
import streams
import argparse
import strutils
import sets
import json
import stats
import math
import bpbiopkg/pedfile
import ./results_html

type Stat4 = object
    dp: RunningStat
    gtdp: RunningStat # depth of genotyped sites
    un: RunningStat
    ab: RunningStat

type count* = object
  nref*: uint32
  nalt*: uint32
  nother*: uint32

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
  #xbounds: array[2, uint16]
  #ybounds: array[2, uint16]
  sites_tested: int
  ibs: seq[uint16]
  n: seq[uint16]
  hets: seq[uint16]
  homs: seq[uint16]
  shared_hom_alts: seq[uint16]
  samples: seq[string]
  # n-samples * n_sites
  allele_counts: seq[seq[count]]
  x_allele_counts: seq[seq[count]]
  y_allele_counts: seq[seq[count]]

proc `%`*(v:uint16): JsonNode =
  new(result)
  result.kind = JInt
  result.num = v.int64

type pair = tuple[a:string, b:string, rel:float64]
proc `%`*(p:pair): JsonNode =
  return %*{"a":p.a, "b":p.b, "rel":p.rel}

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


proc n_samples(r: relation_matrices): int {.inline.} =
  return r.samples.len

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
      if grelatedness > 0.125:
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


type counts* = object
  sample_name*: string
  xbounds*: array[2, uint16]
  ybounds*: array[2, uint16]
  sites*: seq[count]

proc read_extracted(paths: seq[string]): relation_matrices =
  var n_samples = paths.len

  # aggregated from all samples
  result = relation_matrices(ibs: newSeq[uint16](n_samples * n_samples),
                             n: newSeq[uint16](n_samples * n_samples),
                             shared_hom_alts: newSeq[uint16](n_samples * n_samples),
                             hets: newSeq[uint16](n_samples),
                             homs: newSeq[uint16](n_samples),
                             samples: newSeq[string](n_samples),
                             allele_counts: newSeq[seq[count]](n_samples),
                             x_allele_counts: newSeq[seq[count]](n_samples),
                             y_allele_counts: newSeq[seq[count]](n_samples),
                             )
  var
    nsites = 0'u16
    nxsites = 0'u16
    nysites = 0'u16
    last_nsites = 0'u16
    last_nxsites = 0'u16
    last_nysites = 0'u16

  for i, p in paths:
    var f = newFileStream(p, fmRead)
    var sl: uint8 = 0
    discard f.readData(sl.addr, sizeof(sl))
    result.samples[i] = newString(sl)
    discard f.readData(result.samples[i][0].addr, sl.int)
    discard f.readData(nsites.addr, nsites.sizeof.int)
    discard f.readData(nxsites.addr, nxsites.sizeof.int)
    discard f.readData(nysites.addr, nysites.sizeof.int)
    if i > 0:
      doAssert nsites == last_nsites
      doAssert nxsites == last_nxsites
      doAssert nysites == last_nysites

    last_nsites = nsites
    last_nxsites = nxsites
    last_nysites = nysites
    result.allele_counts[i] = newSeq[count](nsites)
    result.x_allele_counts[i] = newSeq[count](nxsites)
    result.y_allele_counts[i] = newSeq[count](nysites)
    doAssert nsites.int * sizeof(result.allele_counts[i][0]) == f.readData(result.allele_counts[i][0].addr, nsites.int * sizeof(result.allele_counts[i][0]))
    doAssert nxsites.int * sizeof(result.x_allele_counts[i][0]) == f.readData(result.x_allele_counts[i][0].addr, nsites.int * sizeof(result.x_allele_counts[i][0]))
    doAssert nysites.int * sizeof(result.y_allele_counts[i][0]) == f.readData(result.y_allele_counts[i][0].addr, nsites.int * sizeof(result.y_allele_counts[i][0]))


    f.close()


proc write(fh:File, sample_names: seq[string], stats: seq[Stat4], gt_counts: array[5, seq[uint16]]) =
  fh.write("#sample\tgt_depth_mean\tgt_depth_sd\tgt_depth_skew\tdepth_mean\tdepth_sd\tdepth_skew\tab_mean\tab_std\tn_hom_ref\tn_het\tn_hom_alt\tn_unknown\tp_middling_ab\n")
  for i, sample in sample_names:
    fh.write(&"{sample}\t")
    fh.write(&"{stats[i].gtdp.mean():.1f}\t{stats[i].gtdp.standard_deviation():.1f}\t{stats[i].gtdp.skewness()}\t")
    fh.write(&"{stats[i].dp.mean():.1f}\t{stats[i].dp.standard_deviation():.1f}\t{stats[i].dp.skewness()}\t")
    fh.write(&"{stats[i].ab.mean():.1f}\t{stats[i].ab.standard_deviation():.1f}\t{gt_counts[0][i]}\t{gt_counts[1][i]}\t{gt_counts[2][i]}\t{gt_counts[3][i]}\t")
    fh.write(&"{gt_counts[4][i].float / (gt_counts[0][i] + gt_counts[1][i] + gt_counts[2][i] + gt_counts[3][i] + gt_counts[4][i]).float:.3f}\n")
  fh.close()

proc toj(samples: seq[string], stats: seq[Stat4], gt_counts: array[5, seq[uint16]]): string =
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
      "n_known": gt_counts[0][i] + gt_counts[1][i] + gt_counts[2][i],
      "p_middling_ab": gt_counts[4][i].float / (gt_counts[0][i] + gt_counts[1][i] + gt_counts[2][i] + gt_counts[3][i] + gt_counts[4][i]).float,
    }
    ))
  result.add("]")

proc rel_main*() =
  ## need to track samples names from bams first, then vcfs since
  ## thats the order for the alts array.
  var argv = commandLineParams()
  if argv[0] == "relate": argv = argv[1..argv.high]

  var p = newParser("somalier relate"):
    help("calculate relatedness among samples from extracted, genotype-like information")
    option("-g", "--groups", help="""optional path  to expected groups of samples (e.g. tumor normal pairs).
specified as comma-separated groups per line e.g.:
    normal1,tumor1a,tumor1b
    normal2,tumor2a""")
    option("-p", "--ped", help="optional path to a ped/fam file indicating the expected relationships among samples.")
    option("-d", "--min-depth", default="7", help="only genotype sites with at least this depth.")
    option("-o", "--output-prefix", help="output prefix for results.", default="somalier")
    arg("extracted", nargs= -1, help="$sample.somalier files for each sample.")


  var opts = p.parse(argv)
  if opts.help:
    quit 0
  if opts.extracted.len == 0:
    echo p.help
    quit "[somalier] specify at least 1 extracted somalier file"
  var
    groups: seq[pair]
    samples: seq[Sample]
    min_depth = parseInt(opts.min_depth)

  if not opts.output_prefix.endswith(".") or opts.output_prefix.endswith("/"):
    opts.output_prefix &= '.'

  var
    final = read_extracted(opts.extracted)
    sample_names = final.samples
    n_samples = sample_names.len

  if opts.ped != "":
    samples = parse_ped(opts.ped)

  stderr.write_line &"[somalier] collected sites from all {n_samples} samples"
  groups.add_ped_samples(samples, final.samples)
  groups.add(readGroups(opts.groups))

  var t0 = cpuTime()
  var nsites = final.allele_counts[0].len
  var n_used_sites = 0
  var alts = newSeq[int8](n_samples)
  # counts of hom-ref, het, hom-alt, unk, hets outside of 0.2..0.8
  var gt_counts : array[5, seq[uint16]]
  var stats = newSeq[Stat4](n_samples)
  for i in 0..<gt_counts.len:
    gt_counts[i] = newSeq[uint16](n_samples)


  #[
  var
    a = newSeq[float32](final.allele_counts[0].len)
    b = newSeq[float32](final.allele_counts[0].len)
  for i in 0..<final.allele_counts.len:
    for k, p in final.allele_counts[i]:
      a[k] = p.ab(min_depth)

    for o in 0..<final.allele_counts.len:
      if i == o: continue
      for k, p in final.allele_counts[i]:
        b[k] = p.ab(min_depth)

      # estimate contamination of a, by b
      var res = estimate_contamination(a, b)
      if res[0] == 0 or res[1] < 10: continue
      echo sample_names[i], " contaminated by ", sample_names[o], " =>", res
  ]#
  for rowi in 0..<nsites:
    var nun = 0
    for i, stat in stats.mpairs:
      var c = final.allele_counts[i][rowi]
      var abi = c.ab(min_depth)
      stat.dp.push(int(c.nref + c.nalt))
      if c.nref > 0'u32 or c.nalt > 0'u32 or c.nother > 0'u32:
        stat.un.push(c.nother.float64 / float64(c.nref + c.nalt + c.nother))
      if c.nref.float > min_depth / 2 and c.nalt.float > min_depth / 2:
        stat.ab.push(abi)
      if abi != -1:
        stat.gtdp.push(int(c.nref + c.nalt))

      alts[i] = abi.alts
      if abi > 0.02 and abi < 0.98 and (abi < 0.2 or abi > 0.8):
        gt_counts[4][i].inc
      if alts[i] == -1:
        nun.inc
        gt_counts[3][i].inc
      else:
        gt_counts[alts[i].int][i].inc

    if nun.float64 / n_samples.float64 > 0.6: continue
    n_used_sites += 1

    discard krelated(alts, final.ibs, final.n, final.hets, final.homs, final.shared_hom_alts, n_samples)

  stderr.write_line &"time to calculate relatedness on {n_used_sites} usable sites: {cpuTime() - t0:.3f}"
  var
    fh_tsv:File
    fh_samples:File
    fh_html:File
    grouped: seq[pair]

  # empty this so it doesn't get sent to html
  final.allele_counts.setLen(0)

  if not open(fh_tsv, opts.output_prefix & "pairs.tsv", fmWrite):
    quit "couldn't open output file"
  if not open(fh_samples, opts.output_prefix & "samples.tsv", fmWrite):
    quit "couldn't open output file"
  if not open(fh_html, opts.output_prefix & "html", fmWrite):
    quit "couldn't open html output file"

  fh_tsv.write_line '#', header.replace("$", "")

  var j = % final
  j["expected-relatedness"] = %* groups

  fh_html.write(tmpl_html.replace("<INPUT_JSON>", $j).replace("<SAMPLE_JSON>", toj(sample_names, stats, gt_counts)))
  fh_html.close()
  stderr.write_line("[somalier] wrote interactive HTML output to: ",  opts.output_prefix & "html")

  fh_samples.write(sample_names, stats, gt_counts)

  for rel in relatedness(final, grouped):
    fh_tsv.write_line $rel

  fh_tsv.close()
  grouped.write(opts.output_prefix)

  stderr.write_line("[somalier] wrote groups to: ",  opts.output_prefix & "groups.tsv")
