import os
import strformat
import random
import bitset
import times
import streams
import algorithm
import argparse
import strutils
import tables
import sets
import json
import ./litestats
import math
import slivarpkg/pedfile
import ./results_html
import ./estimate_contamination

type Stat4 = object
    dp: RunningStat
    gtdp: RunningStat # depth of genotyped sites
    un: RunningStat
    ab: RunningStat

    x_dp: RunningStat
    x_hom_ref: int
    x_het: int
    x_hom_alt: int

    y_dp: RunningStat

type allele_count* = object
  nref*: uint32
  nalt*: uint32
  nother*: uint32

const formatVersion* = 2'u8

type counts* = object
  sample_name*: string
  sites*: seq[allele_count]
  x_sites*: seq[allele_count]
  y_sites*: seq[allele_count]

type relation_matrices = object
  sites_tested: int
  ibs: seq[uint16]
  n: seq[uint16]
  x: seq[uint16]
  shared_hom_alts: seq[uint16]
  samples: seq[string]
  # n-samples * n_sites
  allele_counts: seq[seq[allele_count]]
  x_allele_counts: seq[seq[allele_count]]
  y_allele_counts: seq[seq[allele_count]]
  stats: seq[Stat4]
  gt_counts: array[5, seq[uint16]]
  genotypes: seq[genotypes]
  x_genotypes: seq[genotypes]

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
  x_ibs0: uint16
  x_ibs2: uint16
  n: uint16


## structure for fast small json for plotly
## we use a seq of these separated by the relatedness
## as that's what plotly uses.
type relations = object
  # all sample-pairs with the same expected relatedness are stored together
  # this makes plotly faster
  expected_relatedness: float64
  # TODO: might want to use sample indexes here...
  text: seq[string]
  ibs0: seq[uint16]
  ibs2: seq[uint16]
  shared_hets: seq[uint16]
  shared_hom_alts: seq[uint16]
  concordance: seq[float32]
  relatedness: seq[float32]

proc hom_alt_concordance(r: relation): float64 {.inline.} =
  return (r.shared_hom_alts.float64 - 2 * r.ibs0.float64) / min(r.hom_alts_a, r.hom_alts_b).float64

proc rel(r:relation): float64 {.inline.} =
  return (r.shared_hets.float64 - 2 * r.ibs0.float64) / min(r.hets_a, r.hets_b).float64

proc add*(rt:var seq[relations], rel:relation, expected_relatedness:float) =

  # this keeps order so that unrelateds are first.
  var added: bool = false
  var i = 0
  for r in rt.mitems:
    if abs(r.expected_relatedness - expected_relatedness) < 0.001:
      r.text.add(rel.sample_a & "<br>" & rel.sample_b)
      r.ibs0.add(rel.ibs0)
      r.ibs2.add(rel.ibs2)
      r.shared_hets.add(rel.shared_hets)
      r.shared_hom_alts.add(rel.shared_hom_alts)
      r.concordance.add(rel.hom_alt_concordance)
      r.relatedness.add(rel.rel)
      added = true
      break
    if r.expected_relatedness > expected_relatedness: break
    i += 1

  if not added:
    rt.insert(relations(expected_relatedness: expected_relatedness), i)
    # recurse and add not that we have the correct position.
    rt.add(rel, expected_relatedness)


const header = "$sample_a\t$sample_b\t$relatedness\t$hom_concordance\t$hets_a\t$hets_b\t$shared_hets\t$hom_alts_a\t$hom_alts_b\t$shared_hom_alts\t$ibs0\t$ibs2\t$n\t$x_ibs0\t$x_ibs2\t$expected_relatedness"

proc tsv(r:relation, expected_relatedness:float= -1.0): string =
  result = &"{r.sample_a}\t{r.sample_b}\t{r.rel:.3f}\t{r.hom_alt_concordance:.3f}\t{r.hets_a}\t{r.hets_b}\t{r.shared_hets}\t{r.hom_alts_a}\t{r.hom_alts_b}\t{r.shared_hom_alts}\t{r.ibs0}\t{r.ibs2}\t{r.n}\t{r.x_ibs0}\t{r.xibs2}\t{expected_relatedness}"


proc `%`*(v:uint16): JsonNode =
  result = JsonNode(kind: JInt, num:v.int64)

type pair = tuple[a:string, b:string, rel:float64]
proc `%`*(p:pair): JsonNode =
  return %*{"a":p.a, "b":p.b, "rel":p.rel}

proc cmp_pair(a: pair, b:pair): int =
  result = cmp(a.a, b.a)
  if result == 0: result = cmp(a.b, b.b)

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
  var ss = initHashSet[string]()
  for s in sample_names: ss.incl(s) # use a set for better lookup.
  for i, sampleA in samples[0..<samples.high]:
    if sampleA.id notin ss: continue
    for j, sampleB in samples[i + 1..samples.high]:
      if sampleB.id notin ss: continue
      var rel = sampleA.relatedness(sampleB)
      if rel <= 0: continue
      if sampleA.id < sampleB.id:
        grouped.add((sampleA.id, sampleB.id, rel))
      else:
        grouped.add((sampleB.id, sampleA.id, rel))


proc readGroups(path:string, existing_groups: var seq[pair]): seq[pair] =
  result = newSeq[pair]()
  if path == "":
    return

  var extbl = newTable[string, seq[pair]]()
  # seen makes sure we don't add a pair that's already present
  var seen = newTable[tuple[a:string, b:string], bool]()
  for ex in existing_groups:
    extbl.mgetOrPut(ex.a, @[]).add(ex)
    extbl.mgetOrPut(ex.b, @[]).add(ex)
    seen[(ex.a, ex.b)] = true

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

        var added = result[result.high]
        for up in extbl.getOrDefault(added.a, @[]):
          var toadd:pair = (added.b, up.b, up.rel)
          if toadd.b < toadd.a: swap(toadd.a, toadd.b)
          # we know up.a and added.a ar already pairs so we need to pair up.b
          # and added.a
          if (toadd.a, toadd.b) notin seen:
            seen[(toadd.a, toadd.b)] = true
            existing_groups.add(toadd)

        # repeat above for b
        for up in extbl.getOrDefault(added.b, @[]):
          var toadd:pair = (added.a, up.b, up.rel)
          if toadd.b < toadd.a: swap(toadd.a, toadd.b)
          if (toadd.a, toadd.b) notin seen:
            seen[(toadd.a, toadd.b)] = true
            existing_groups.add(toadd)


proc n_samples(r: relation_matrices): int {.inline.} =
  return r.samples.len

iterator relatedness(r: var relation_matrices, grouped: var seq[pair]): relation =
  var
    ir:IBSResult
    xir:IBSResult
    sample_names = r.samples

  for j in 0..<r.genotypes.high:
    let hets_j = r.gt_counts[1][j]
    for k in (j + 1) .. r.genotypes.high:
      let hets_k = r.gt_counts[1][k]

      var bottom: float64 = 0
      if (hets_j > 0'u16) and (hets_k > 0'u16):
        bottom = 2 / (1 / hets_k.float64 + 1 / hets_j.float64).float64
      else:
        bottom = max(hets_j, hets_k).float64
      if bottom == 0:
        # can't calculate relatedness
        bottom = -1'f64

      ir = r.genotypes[j].IBS(r.genotypes[k])
      let grel = (ir.shared_hets.float64 - 2 * ir.IBS0.float64) / bottom
      if grel > 0.125:
        grouped.add((sample_names[j], sample_names[k], grel))
      # now fill the matrices so they can be used from javascript
      r.ibs[j * r.n_samples + k] = ir.IBS0.uint16
      r.ibs[k * r.n_samples + j] = ir.shared_hets.uint16
      r.n[j * r.n_samples + k] = ir.N.uint16
      r.n[k * r.n_samples + j] = ir.IBS2.uint16
      r.shared_hom_alts[j * r.n_samples + k] = ir.shared_hom_alts.uint16

      xir = r.x_genotypes[j].XIBS(r.x_genotypes[k])
      r.x[j * r.n_samples + k] = xir.IBS0.uint16
      r.x[k * r.n_samples + j] = xir.IBS2.uint16

      yield relation(sample_a: sample_names[j],
                     sample_b: sample_names[k],
                     hets_a: hets_j, hets_b: hets_k,
                     hom_alts_a: r.gt_counts[2][j], hom_alts_b: r.gt_counts[2][k],
                     ibs0: ir.IBS0.uint16,
                     shared_hets: ir.shared_hets.uint16,
                     shared_hom_alts: ir.shared_hom_alts.uint16,
                     ibs2: ir.IBS2.uint16,
                     n: ir.N.uint16,
                     x_ibs0: xir.IBS0.uint16,
                     x_ibs2: xir.IBS2.uint16,
                     )


template proportion_other(c:allele_count): float =
  if c.nother == 0: 0'f else: c.nother.float / (c.nother + c.nref + c.nalt).float

proc ab*(c:allele_count, min_depth:int): float {.inline.} =
  ## get the allele balance for the allele_count object while requing a min depth
  # allow high-ish proportion other see: 
  # https://github.com/brentp/somalier/issues/26#issuecomment-543120582
  if c.proportion_other > 0.1: return -1
  if int(c.nref + c.nalt) < min_depth:
    return -1
  if c.nalt == 0:
    return 0
  result = c.nalt.float / (c.nalt + c.nref).float

proc alts*(ab:float, ab_cutoff:float=0.04): int8 {.inline.} =
  if ab < 0: return -1
  if ab < ab_cutoff: return 0
  if ab > (1 - ab_cutoff): return 2
  if ab >= 0.2 and ab <= 0.8: return 1
  return -1

{.push checks: off, optimization: speed.}
template depth*(c:allele_count): uint32 =
  c.nref + c.nalt

proc alts*(c:allele_count, min_depth:int=7): int8 {.inline.} =
  if c.proportion_other > 0.04: return -1
  if int(c.nref + c.nalt) < min_depth: return -1
  if c.nref == 0: return 2
  if c.nalt == 0: return 0
  var ab = c.nalt.float / (c.depth).float
  return ab.alts

proc to_bits*(cs:seq[allele_count], min_depth:int=7): genotypes =
  result.hom_ref = create_bitset(cs.len)
  result.het = create_bitset(cs.len)
  result.hom_alt = create_bitset(cs.len)

  for i, c in cs:
    var a = c.alts(min_depth)
    if a == 0:
      result.hom_ref.set(i)
    elif a == 1:
      result.het.set(i)
    elif a == 2:
      result.hom_alt.set(i)
{.pop.}

proc fill_sample_info(r:var relation_matrices, sample_i:int, min_depth:int, unk2hr:bool) =

  var n = r.allele_counts[sample_i].len
  r.genotypes[sample_i].hom_ref = create_bitset(n)
  r.genotypes[sample_i].het = create_bitset(n)
  r.genotypes[sample_i].hom_alt = create_bitset(n)

  var stat = r.stats[sample_i]
  for k, c in r.allele_counts[sample_i]:
    var abi = c.ab(min_depth)
    if abi < 0 and unk2hr: abi = 0
    stat.dp.push(int(c.nref + c.nalt))
    if c.nref > 0'u32 or c.nalt > 0'u32 or c.nother > 0'u32:
      stat.un.push(c.nother.float64 / float64(c.nref + c.nalt + c.nother))
      # TODO: why is this here?
    if c.nref.float > min_depth / 2 or c.nalt.float > min_depth / 2:
      stat.ab.push(abi)
    if abi != -1:
      stat.gtdp.push(int(c.nref + c.nalt))
    var alt = abi.alts
    if alt < 0 and unk2hr: alt = 0
    if abi > 0.02 and abi < 0.98 and (abi < 0.1 or abi > 0.9):
      r.gt_counts[4][sample_i].inc
    if alt == -1:
      r.gt_counts[3][sample_i].inc
    else:
      r.gt_counts[alt][sample_i].inc
      if alt == 0:
        r.genotypes[sample_i].hom_ref.set(k)
      elif alt == 1:
        r.genotypes[sample_i].het.set(k)
      elif alt == 2:
        r.genotypes[sample_i].hom_alt.set(k)


  n = r.x_allele_counts[sample_i].len
  r.x_genotypes[sample_i].hom_ref = create_bitset(n)
  r.x_genotypes[sample_i].het = create_bitset(n)
  r.x_genotypes[sample_i].hom_alt = create_bitset(n)
  for k, c in r.x_allele_counts[sample_i]:
    var alt = c.alts
    if alt == -1: continue
    stat.x_dp.push(c.depth.float)
    if alt == 0:
      stat.x_hom_ref.inc
      r.x_genotypes[sample_i].hom_ref.set(k)
    elif alt == 1:
      stat.x_het.inc
      r.x_genotypes[sample_i].het.set(k)
    elif alt == 2:
      stat.x_hom_alt.inc
      r.x_genotypes[sample_i].hom_alt.set(k)

  for c in r.y_allele_counts[sample_i]:
    var alt = c.alts
    if alt == -1: continue
    stat.y_dp.push(c.depth.float)
  r.stats[sample_i] = stat


proc read_extracted*(paths: seq[string], min_depth:int, unk2hr:bool): relation_matrices =
  var n_samples = paths.len

  # aggregated from all samples
  result = relation_matrices(ibs: newSeq[uint16](n_samples * n_samples),
                             n: newSeq[uint16](n_samples * n_samples),
                             shared_hom_alts: newSeq[uint16](n_samples * n_samples),
                             x: newSeq[uint16](n_samples * n_samples),
                             samples: newSeq[string](n_samples),
                             allele_counts: newSeq[seq[allele_count]](n_samples),
                             x_allele_counts: newSeq[seq[allele_count]](n_samples),
                             y_allele_counts: newSeq[seq[allele_count]](n_samples),
                             stats: newSeq[Stat4](n_samples),
                             genotypes: newSeq[genotypes](n_samples),
                             x_genotypes: newSeq[genotypes](n_samples),
                             )
  var
    nsites = 0'u16
    nxsites = 0'u16
    nysites = 0'u16
    last_nsites = 0'u16
    last_nxsites = 0'u16
    last_nysites = 0'u16

  for i in 0..<result.gt_counts.len:
    result.gt_counts[i] = newSeq[uint16](n_samples)

  for i, p in paths:
    var f = newFileStream(p, fmRead)
    var sl: uint8 = 0
    discard f.readData(sl.addr, sizeof(sl))
    doAssert sl == formatVersion, &"expected matching versions got {sl}, expected {formatVersion}"
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
    result.allele_counts[i] = newSeq[allele_count](nsites)
    result.x_allele_counts[i] = newSeq[allele_count](nxsites)
    result.y_allele_counts[i] = newSeq[allele_count](nysites)
    if nsites > 0'u16:
      doAssert nsites.int * sizeof(result.allele_counts[i][0]) == f.readData(result.allele_counts[i][0].addr, nsites.int * sizeof(result.allele_counts[i][0])), &"error in file: {p}"
    if nxsites > 0'u16:
      doAssert nxsites.int * sizeof(result.x_allele_counts[i][0]) == f.readData(result.x_allele_counts[i][0].addr, nxsites.int * sizeof(result.x_allele_counts[i][0])), &"error in file: {p}"
    if nysites > 0'u16:
      doAssert nysites.int * sizeof(result.y_allele_counts[i][0]) == f.readData(result.y_allele_counts[i][0].addr, nysites.int * sizeof(result.y_allele_counts[i][0])), &"error in file: {p}"

    f.close()

    result.fill_sample_info(i, min_depth, unk2hr)


proc write(fh:File, sample_names: seq[string], stats: seq[Stat4], gt_counts: array[5, seq[uint16]], sample_sex: TableRef[string, string]) =
  fh.write("#sample\tpedigree_sex\tgt_depth_mean\tgt_depth_sd\tdepth_mean\tdepth_sd\tab_mean\tab_std\tn_hom_ref\tn_het\tn_hom_alt\tn_unknown\tp_middling_ab\t")
  fh.write("X_depth_mean\tX_n\tX_hom_ref\tX_het\tX_hom_alt\t")
  fh.write("Y_depth_mean\tY_n\n")
  for i, sample in sample_names:
    fh.write(&"{sample}\t{sample_sex.getOrDefault(sample)}\t")
    fh.write(&"{stats[i].gtdp.mean:.1f}\t{stats[i].gtdp.standard_deviation():.1f}\t")
    fh.write(&"{stats[i].dp.mean:.1f}\t{stats[i].dp.standard_deviation():.1f}\t")
    fh.write(&"{stats[i].ab.mean:.2f}\t{stats[i].ab.standard_deviation():.2f}\t{gt_counts[0][i]}\t{gt_counts[1][i]}\t{gt_counts[2][i]}\t{gt_counts[3][i]}\t")
    fh.write(&"{gt_counts[4][i].float / (gt_counts[0][i] + gt_counts[1][i] + gt_counts[2][i] + gt_counts[3][i] + gt_counts[4][i]).float:.3f}\t")
    fh.write(&"{stats[i].x_dp.mean:.2f}\t{stats[i].x_dp.n}\t{stats[i].x_hom_ref}\t{stats[i].x_het}\t{stats[i].x_hom_alt}\t")
    fh.write(&"{stats[i].y_dp.mean:.2f}\t{stats[i].y_dp.n}\n")
  fh.close()

proc toj(sample_names: seq[string], stats: seq[Stat4], gt_counts: array[5, seq[uint16]], sample_sex: TableRef[string, string]): string =
  result = newStringOfCap(10000)
  result.add("[")
  for i, s in sample_names:
    if i > 0: result.add(",\n")
    result.add($(%* {
      "sample": s,
      "sex": sample_sex.getOrDefault(s, "unknown"),

      "gt_depth_mean": stats[i].gtdp.mean,

      "depth_mean": stats[i].dp.mean,

      "ab_mean": stats[i].ab.mean,
      "pct_other_alleles": 100.0 * stats[i].un.mean,
      "n_hom_ref": gt_counts[0][i],
      "n_het": gt_counts[1][i],
      "n_hom_alt": gt_counts[2][i],
      "n_unknown": gt_counts[3][i],
      "n_known": gt_counts[0][i] + gt_counts[1][i] + gt_counts[2][i],
      "p_middling_ab": gt_counts[4][i].float / (gt_counts[0][i] + gt_counts[1][i] + gt_counts[2][i] + gt_counts[3][i] + gt_counts[4][i]).float,

      "x_depth_mean": 2 * stats[i].x_dp.mean / stats[i].gtdp.mean,
      "x_hom_ref": stats[i].x_hom_ref,
      "x_het": stats[i].x_het,
      "x_hom_alt": stats[i].x_hom_alt,

      "y_depth_mean": 2 * stats[i].y_dp.mean / stats[i].gtdp.mean,
    }
    ))
  result.add("]")

proc to_sex_lookup(samples: seq[Sample]): TableRef[string, string] =
  result = newTable[string, string]()
  for s in samples:
    result[s.id] = if s.sex == 1: "male" elif s.sex == 2: "female" else: "unknown"

proc rel_main*() =
  ## need to track samples names from bams first, then vcfs since
  ## thats the order for the alts array.
  randomize()
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
    flag("-u", "--unknown", help="set unknown genotypes to hom-ref. it is often preferable to use this with VCF samples that were not jointly called")
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
    unk2hr = opts.unknown

  if not opts.output_prefix.endswith(".") or opts.output_prefix.endswith("/"):
    opts.output_prefix &= '.'

  var t0 = cpuTime()
  var final = read_extracted(opts.extracted, min_depth, unk2hr)
  var n_samples = final.samples.len
  stderr.write_line &"[somalier] time to read files and get per-sample stats for {n_samples} samples: {cpuTime() - t0:.2f}"
  t0 = cpuTime()

  if opts.ped != "":
    samples = parse_ped(opts.ped)
  if samples.len > 30_000:
    stderr.write_line "[somalier] WARNING!! somalier will work fine for even 100K samples, but it is not optimal for such scenarios."
    stderr.write_line "[somalier] ......... please open an issue at: https://github.com/brentp/somalier/issues as larger cohorts"
    stderr.write_line "[somalier] ......... can be supported."

  groups.add_ped_samples(samples, final.samples)
  # send in groups so we can adjust baed on self-self samples
  groups.add(readGroups(opts.groups, groups))
  stderr.write_line &"[somalier] time to get expected relatedness from pedigree graph: {cpuTime() - t0:.2f}"


  var
    fh_tsv:File
    fh_samples:File
    fh_html:File
    grouped: seq[pair]

  # empty this so it doesn't get sent to html

  if not open(fh_tsv, opts.output_prefix & "pairs.tsv", fmWrite):
    quit "couldn't open output file"
  if not open(fh_samples, opts.output_prefix & "samples.tsv", fmWrite):
    quit "couldn't open output file"
  if not open(fh_html, opts.output_prefix & "html", fmWrite):
    quit "couldn't open html output file"

  t0 = cpuTime()

  var tmpls = tmpl_html.split("<INPUT_JSON>")
  var sample_sexes = samples.to_sex_lookup
  fh_html.write(tmpls[0].replace("<SAMPLE_JSON>", toj(final.samples, final.stats, final.gt_counts, sample_sexes)))

  var rels: seq[relations]

  var proportion_sampled = 400_000'f64 / float64(final.samples.len * max(1, final.samples.len - 1))
  if proportion_sampled < 1:
    stderr.write_line &"[somalier] html output will have unrelated sample-pairs subset to {100 * proportion_sampled:.2f}% of points"

  fh_tsv.write_line '#', header.replace("$", "")
  var npairs:int
  var nrels:int
  sort(groups, cmp_pair)
  for rel in final.relatedness(grouped):
    var idx = groups.binarySearch((rel.sample_a, rel.sample_b, -1.0), cmp_pair)
    if idx == -1:
      idx = groups.binarySearch((rel.sample_b, rel.sample_a, -1.0), cmp_pair)
    let expected_relatedness = if idx == -1: -1'f else: groups[idx].rel
    if (expected_relatedness != -1) or (rand(1'f32) < proportion_sampled) or rel.rel > 0.08:
      rels.add(rel, max(0, expected_relatedness))
      nrels.inc

    fh_tsv.write_line rel.tsv(expected_relatedness)
    npairs.inc
  stderr.write_line &"[somalier] time to calculate all vs all relatedness for all {npairs} combinations: {cpuTime() - t0:.2f}"

  fh_html.write(%* rels)
  fh_html.write(tmpls[1])
  fh_html.close()
  stderr.write_line(&"[somalier] wrote interactive HTML output for {nrels} pairs to: ",  opts.output_prefix & "html")

  fh_samples.write(final.samples, final.stats, final.gt_counts, sample_sexes)

  fh_tsv.close()
  grouped.write(opts.output_prefix)

  stderr.write_line("[somalier] wrote groups to: ",  opts.output_prefix & "groups.tsv")
  stderr.write_line("[somalier] wrote samples to: ",  opts.output_prefix & "samples.tsv")
  stderr.write_line("[somalier] wrote pair-wise relatedness metrics to: ",  opts.output_prefix & "pairs.tsv")
