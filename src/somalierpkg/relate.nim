import os
import strformat
import random
import bitset
import times
import streams
import algorithm
import argparse
import sequtils
import strutils
import tables
import sets
import json
import ./litestats
import math
import pedfile
import ./results_html
import ./common
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


type relation_matrices = object {.shallow.}
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
  het_ab: uint16
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
  n: seq[uint16]

proc hom_alt_concordance(r: relation): float64 {.inline.} =
  return (r.shared_hom_alts.float64 - 2 * r.ibs0.float64) / max(1'u16, min(r.hom_alts_a, r.hom_alts_b)).float64

proc rel(r:relation): float64 {.inline.} =
  return 2 * (r.shared_hets.float64 - 2 * r.ibs0.float64) / r.het_ab.float64

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
      r.n.add(rel.n)
      added = true
      break
    if r.expected_relatedness > expected_relatedness: break
    i += 1

  if not added:
    rt.insert(relations(expected_relatedness: expected_relatedness), i)
    # recurse and add now that we have the correct position.
    rt.add(rel, expected_relatedness)


const header = "$sample_a\t$sample_b\t$relatedness\t$ibs0\t$ibs2\t$hom_concordance\t$hets_a\t$hets_b\t$hets_ab\t$shared_hets\t$hom_alts_a\t$hom_alts_b\t$shared_hom_alts\t$n\t$x_ibs0\t$x_ibs2\t$expected_relatedness"

proc tsv(r:relation, expected_relatedness:float= -1.0): string =
  result = &"{r.sample_a}\t{r.sample_b}\t{r.rel:.3f}\t{r.ibs0}\t{r.ibs2}\t{r.hom_alt_concordance:.3f}\t{r.hets_a}\t{r.hets_b}\t{r.het_ab}\t{r.shared_hets}\t{r.hom_alts_a}\t{r.hom_alts_b}\t{r.shared_hom_alts}\t{r.n}\t{r.x_ibs0}\t{r.xibs2}\t{expected_relatedness}"


proc to_sex_lookup(samples: seq[Sample]): TableRef[string, string] =
  result = newTable[string, string]()
  for s in samples:
    result[s.id] = if s.sex == 1: "male" elif s.sex == 2: "female" else: "unknown"


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
      if grp.rel > 0.98:
        fh_groups.write(&"{grp.a},{grp.b}\t1\n")
      else:
        fh_groups.write(&"{grp.a},{grp.b}\t{grp.rel:.2f}\n")
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
      rel = round(parseFloat(tmp[1]), 2)
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

proc relatedness(r: var relation_matrices, j: int, k:int): relation {.inline.} =
  var j = j
  var k = k
  if j > k:
    let tmp = j
    j = k
    k = tmp
  let hets_k = r.gt_counts[1][k]
  let hets_j = r.gt_counts[1][j]

  if r.n[j * r.n_samples + k] > 0'u16: # used previously calculated data
    return relation(#sample_a: sample_names[j],
                 #sample_b: sample_names[k],
                 hets_a: hets_j, hets_b: hets_k,
                 hom_alts_a: r.gt_counts[2][j], hom_alts_b: r.gt_counts[2][k],
                 ibs0: r.ibs[j * r.n_samples + k],
                 shared_hets: r.ibs[k * r.n_samples + j],
                 shared_hom_alts: r.shared_hom_alts[j * r.n_samples + k],
                 het_ab: r.shared_hom_alts[k * r.n_samples + j],
                 ibs2: r.n[k * r.n_samples + j],
                 n: r.n[j * r.n_samples + k],
                 x_ibs0: r.x[j * r.n_samples + k],
                 x_ibs2: r.x[k * r.n_samples + j],
               )

  let ir = r.genotypes[j].IBS(r.genotypes[k])
  # now fill the matrices so they can be used from javascript
  r.ibs[j * r.n_samples + k] = ir.IBS0.uint16
  r.ibs[k * r.n_samples + j] = ir.shared_hets.uint16
  r.n[j * r.n_samples + k] = ir.N.uint16
  r.n[k * r.n_samples + j] = ir.IBS2.uint16
  r.shared_hom_alts[j * r.n_samples + k] = ir.shared_hom_alts.uint16
  r.shared_hom_alts[k * r.n_samples + j] = ir.het_ab.uint16

  let xir = r.x_genotypes[j].XIBS(r.x_genotypes[k])
  r.x[j * r.n_samples + k] = xir.IBS0.uint16
  r.x[k * r.n_samples + j] = xir.IBS2.uint16

  result = relation(#sample_a: sample_names[j],
                 #sample_b: sample_names[k],
                 hets_a: hets_j, hets_b: hets_k,
                 hom_alts_a: r.gt_counts[2][j], hom_alts_b: r.gt_counts[2][k],
                 ibs0: ir.IBS0.uint16,
                 shared_hets: ir.shared_hets.uint16,
                 shared_hom_alts: ir.shared_hom_alts.uint16,
                 ibs2: ir.IBS2.uint16,
                 n: ir.N.uint16,
                 het_ab: ir.het_ab.uint16,
                 x_ibs0: xir.IBS0.uint16,
                 x_ibs2: xir.IBS2.uint16,
                 )

iterator relatedness(r: var relation_matrices, grouped: var seq[pair]): tuple[r:relation, i:int, j:int] =

  let sample_names = r.samples
  for j in 0..<r.genotypes.high:
    let sample_a = sample_names[j]
    for k in (j + 1) .. r.genotypes.high:
      var r = r.relatedness(j, k)
      r.sample_a = sample_a
      r.sample_b = sample_names[k]
      if r.rel > 0.125:
        grouped.add((r.sample_a, r.sample_b, r.rel))
      yield (r, j, k)


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

var ab_cutoff: float = 0.01
try:
  ab_cutoff = parseFloat(getEnv("SOMALIER_AB_HOM_CUTOFF"))
  if ab_cutoff > 0.5:
    stderr.writeline("[somalier] error setting SOMALIER_AB_HOM_CUTOFF to:" & getEnv("SOMALIER_AB_HOM_CUTOFF"))
    ab_cutoff = 0.01
except:
  discard

proc alts*(ab:float, min_ab:float, ab_cutoff:float=ab_cutoff): int8 {.inline.} =
  if ab < 0: return -1
  if ab < ab_cutoff: return 0
  if ab > (1 - ab_cutoff): return 2
  if ab >= min_ab and ab <= (1 - min_ab): return 1
  return -1

{.push checks: off, optimization: speed.}
template depth*(c:allele_count): uint32 =
  c.nref + c.nalt

proc alts*(c:allele_count, min_ab:float, min_depth:int=7): int8 {.inline.} =
  if c.proportion_other > 0.04: return -1
  if int(c.nref + c.nalt) < min_depth: return -1
  if c.nref == 0: return 2
  if c.nalt == 0: return 0
  var ab = c.nalt.float / (c.depth).float
  return ab.alts(min_ab)

{.pop.}

proc fill_sample_info(r:var relation_matrices, sample_i:int, min_ab:float, min_depth:int, unk2hr:bool) =

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
    var alt = abi.alts(min_ab)
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
    var alt = c.alts(min_ab)
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
    var alt = c.alts(min_ab)
    if alt == -1: continue
    stat.y_dp.push(c.depth.float)
  r.stats[sample_i] = stat


proc read_extracted*(paths: seq[string], min_ab:float, min_depth:int, unk2hr:bool): relation_matrices =
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
    if f == nil:
      raise newException(IOError, "could not open file:" & p)
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

    result.fill_sample_info(i, min_ab, min_depth, unk2hr)

const missing = [".", "0", "-9", ""]

proc high_quality(gt_counts: array[5, seq[uint16]], i:int): bool {.inline.} =
  # less than 3% of sites with allele balance outside of 0.1 .. 0.9
  result = gt_counts[4][i].float / (gt_counts[0][i] + gt_counts[1][i] + gt_counts[2][i] + gt_counts[3][i] + gt_counts[4][i]).float < 0.06
  if not result:
    return false
  result = gt_counts[2][i].float / gt_counts[1][i].float > 0.7

proc samples_have_y_depth(stats: seq[Stat4]): bool =
  var n = 0
  for s in stats:
    n += int(s.y_dp.n > 0)
  return n > 5 or n.float / stats.len.float > 0.1

type SampleLooker = object
  sample_names: seq[string]
  sample_sex: TableRef[string, string]
  sample_table: TableRef[string, Sample]
  pairs: TableRef[string, seq[string]]
  sib_pairs: TableRef[string, seq[string]]
  has_y: bool
  changed_samples: HashSet[string]


proc unrelated(final: var relation_matrices, L: SampleLooker, possible_parents: var seq[string], stats: seq[Stat4], check_parent_sexes:bool=true, level:float=0.06): bool =
  doAssert possible_parents.len == 2, "[somalier] ERROR expected only 2 parents in call to 'unrelated'"

  var
    i = L.sample_table[possible_parents[0]].i
    j = L.sample_table[possible_parents[1]].i

  if i < 0 or j < 0: return false


  let rel = final.relatedness(i, j)
  if rel.rel > level: return false

  # now order possible_parents so that father is first.
  if check_parent_sexes:
    var imale = stats[i].x_het / stats[i].x_hom_alt < 0.05
    var jmale = stats[j].x_het / stats[j].x_hom_alt < 0.05
    if imale == jmale: return false

    if jmale:
      let tmp = possible_parents[0]
      possible_parents[0] = possible_parents[1]
      possible_parents[1] = tmp

  return true

proc related(final: var relation_matrices, L: SampleLooker, possible_parents: var seq[string], stats: seq[Stat4], check_parent_sexes:bool=true, level:float=0.36): bool =
  return not unrelated(final, L, possible_parents, stats, check_parent_sexes, level)


proc add_parents_and_check_sex(final:var relation_matrices, stats: seq[Stat4], gt_counts: array[5, seq[uint16]], i:int, L:var SampleLooker) =
  # first pass here updates sample parents as needed.
  if not gt_counts.high_quality(i):
    return
  let sample_name = L.sample_names[i]
  var sample = L.sample_table.getOrDefault(sample_name, Sample(id: sample_name, family_id: sample_name, sex: -9, phenotype:"-9", maternal_id:"-9", paternal_id:"-9"))
  ## can re-set values for high quality samples based on info in VCF
  if stats[i].x_het / stats[i].x_hom_alt < 0.05 and stats[i].x_dp.n > 10:
    if sample.sex != 1:
      if sample.sex == 2:
        stderr.write_line &"[somalier] setting sex to male for {sample.id}"
      sample.sex = 1
  elif stats[i].x_het / stats[i].x_hom_alt > 0.4 and stats[i].x_dp.n > 10:
    if sample.sex != 2:
      if sample.sex == 1:
        stderr.write_line &"[somalier] setting sex to female for {sample.id}"
      sample.sex = 2
  if L.has_y and sample.sex == 1 and 2 * stats[i].y_dp.mean / stats[i].gtdp.mean < 0.4:
    stderr.write_line &"[somalier] NOTE: apparent loss of Y for {sample.id} with low het-ratio on X chromosome"
  if L.has_y and sample.sex == 2 and 2 * stats[i].y_dp.mean / stats[i].gtdp.mean > 0.4:
    stderr.write_line &"[somalier] NOTE: apparent Y for {sample.id} with high het-ratio on X chromosome"
    sample.sex = -2
  # now look up in parent child pairs. if there are 2 samples, we check
  # that those 2 samples are unrelated. if so, they are mom and dad.
  var possible_parents = L.pairs.getOrDefault(sample.id, @[])
  # TODO: for 3 gens, often won't have exactly 2. need to find exactly 2 that
  # are unrelated...
  # call to unrelated also orders parents so that order is dad, mom as in
  # pedigree file.
  if possible_parents.len == 2 and final.unrelated(L, possible_parents, stats):
    if sample.paternal_id != possible_parents[0]:
      stderr.write_line &"[somalier] NOTE: updating paternal_id for {sample.id} to {possible_parents[0]}"
      sample.paternal_id = possible_parents[0]
    if sample.maternal_id != possible_parents[1]:
      stderr.write_line &"[somalier] NOTE: updating maternal_id for {sample.id} to {possible_parents[1]}"
      sample.maternal_id = possible_parents[1]
  L.sample_table[sample.id] = sample

proc add_parent_to_sibs(final:var relation_matrices, stats: seq[Stat4], gt_counts: array[5, seq[uint16]], L:var SampleLooker) =
  # now, if we have multiple siblings that all have low IBS0 to an additional
  # $sample, we can assume $sample is the (single) parent
  var byParent = newTable[string, seq[string]]()

  for sid, sample in L.sample_table:
    if sample.paternal_id notin missing:
      byParent.mgetOrPut(sample.paternal_id, @[]).add(sample.id)
    if sample.maternal_id notin missing:
      byParent.mgetOrPut(sample.maternal_id, @[]).add(sample.id)

  for parent_id, kid_ids in byParent:
    if kid_ids.len < 2: continue
    var parent_ids = initCountTable[string]()
    for i, k in kid_ids:
      if k notin L.pairs:
        break
      for p in L.pairs[k]:
        parent_ids.inc(p)
    # now, there there should be 1 entry in parent_ids with count == kid_ids.len
    for parent_id, c in parent_ids:
      if c == kid_ids.len:
        let parent = L.sample_table.getOrDefault(parent_id, Sample())
        if parent.id == "": continue
        for kid_id in kid_ids:
          let kid = L.sample_table[kid_id]
          if parent.sex == 1 and kid.paternal_id != parent.id:
              if not kid.paternal_id.endswith("_somalier"):
                stderr.write_line &"[somalier] WARNING: updating paternal id for sample {kid.id} from {kid.paternal_id} to {parent.id}"
              kid.paternal_id = parent.id
          elif parent.sex == 2 and kid.maternal_id != parent.id:
              if not kid.maternal_id.endswith("_somalier"):
                stderr.write_line &"[somalier] WARNING: updating maternal id for sample {kid.id} from {kid.maternal_id} to {parent.id}"
              kid.maternal_id = parent.id

proc remove_spurious_parent_ids(final:var relation_matrices, L:SampleLooker, stats: seq[Stat4]) =
  for id, sample in L.sample_table.mpairs:
    if sample.id in L.sample_table and sample.paternal_id in L.sample_table and not sample.paternal_id.endswith("_somalier"):
      var pair = @[sample.id, sample.paternal_id]
      if not final.related(L, pair, stats, false, 0.33):
        stderr.write_line &"[somalier] removing assigned father from {sample.id} and setting to unknown"
        sample.paternal_id = "-9"
    if sample.id in L.sample_table and sample.maternal_id in L.sample_table and not sample.maternal_id.endswith("_somalier"):
      var pair = @[sample.id, sample.maternal_id]
      if not final.related(L, pair, stats, false, 0.33):
        stderr.write_line &"[somalier] removing assigned mother from {sample.id} and setting to unknown"
        sample.maternal_id = "-9"
    L.sample_table[id] = sample


proc add_siblings(final:var relation_matrices, stats: seq[Stat4], gt_counts: array[5, seq[uint16]], L:var SampleLooker) =
  for sample_name, sib_names in L.sib_pairs:
    let isample = L.sample_table[sample_name]
    let iset = sib_names.toSet
    let i = isample.i
    if i >= 0 and not gt_counts.high_quality(i): continue
    var ipids = [isample.paternal_id, isample.maternal_id]
    let parent_order = ["dad", "mom"]
    for sn in sib_names:

      let jset = L.sib_pairs.getOrDefault(sn, @[]).toSet
      # require that they share the same siblings
      if iset.symmetricDifference(jset).len != 2:
        continue

      var j: int
      var jsample:Sample
      try:
        jsample = L.sample_table[sn]
        j = jsample.i
      except KeyError:
        continue
      if j >= 0 and not gt_counts.high_quality(j): continue
      #if isample.paternal_id notin missing and isample.maternal_id notin missing and isample.paternal_id == jsample.paternal_id and isample.maternal_id == jsample.maternal_id: continue
      ## TODO: some logic problems below. check in CEPH
      var jpids = [jsample.paternal_id, jsample.maternal_id]

      for k, ipid in ipids.mpairs:
        if ipid in [isample.id, jsample.id]:
          continue
        var jpid = jpids[k]

        if jpid in [isample.id, jsample.id]:
          continue

        var changed = false

        if ipid notin missing and jpid in missing:
          jpids[k] = ipid
          jpid = ipid
          changed = true
        elif jpid notin missing and ipid in missing:
          ipid = jpid
          ipids[k] = jpid
          changed = true

        elif ipid in missing and jpid in missing:
          # make a fake dad
          ipid = &"""{parent_order[k]}_{isample.family_id}_somalier"""
          jpid = ipid
          ipids[k] = ipid
          jpids[k] = ipid
          changed = true

        elif ipid != jpid: # both samples had a different parent specified
          if jsample.id notin missing or (final.relatedness(jsample.i, isample.i).rel > 0.42 and final.relatedness(jsample.i, L.sample_table[ipid].i).ibs0 < final.relatedness(jsample.i, L.sample_table[jpid].i).ibs0):
            stderr.write_line &"[somalier] NOTE: apparent siblings {jsample.id} and {isample.id} have a different {parent_order[k]} setting to {ipid} ({jsample.id} had {jpid})"
            jpid = ipid
            jpids[k] = ipid
            changed = true

        if changed:
          if isample.family_id != jsample.family_id:
            stderr.write_line &"[somalier] NOTE updating family id of sample {jsample.id} to sibling {isample.family_id}"
            jsample.family_id = isample.family_id

        if ipid != jpid:
          if final.relatedness(jsample.i, isample.i).rel > 0.42:
            stderr.write_line &"[somalier] ERROR not specified as sibs:", isample, " ", jsample

      if jpids[0] != jsample.id:
        if jpids[0] != jsample.paternal_id:
          var pair = @[jpids[0], jsample.id]
          if pair[0] notin L.sample_table or pair[1] notin L.sample_table or final.related(L, pair, stats, false, level=0.4):
            jsample.paternal_id = jpids[0]

      if jpids[1] != jsample.id:
        if jpids[1] != jsample.maternal_id:
          # make sure we don't accidently set wife as mother just because they
          # share an offspring.
          var pair = @[jpids[1], jsample.id]
          if pair[0] notin L.sample_table or pair[1] notin L.sample_table or final.related(L, pair, stats, false, level=0.4):
            jsample.maternal_id = jpids[1]
      L.sample_table[sn] = jsample

      if ipids[0] != isample.id:
        if ipids[0] != isample.paternal_id:
          var pair = @[ipids[0], isample.id]
          if pair[0] notin L.sample_table or pair[1] notin L.sample_table or final.related(L, pair, stats, false, level=0.4):
            isample.paternal_id = ipids[0]

      if ipids[1] != isample.id:
        if ipids[1] != isample.maternal_id:
          var pair = @[ipids[1], isample.id]
          if pair[0] notin L.sample_table or pair[1] notin L.sample_table or final.related(L, pair, stats, false, level=0.4):
            isample.maternal_id = ipids[1]

      L.sample_table[sample_name] = isample


proc update_family_ids(final:relation_matrices, stats: seq[Stat4], gt_counts: array[5, seq[uint16]], i:int, L:var SampleLooker) =
  # 2nd pass here updates family ids
  let sample_name = L.sample_names[i]
  if sample_name notin L.sample_table:
    let s = Sample(id: sample_name, family_id: sample_name, sex: -9, phenotype:"-9", maternal_id:"-9", paternal_id:"-9")
  var sample = L.sample_table[sample_name]
  if not gt_counts.high_quality(i): return

  let pids = [sample.paternal_id, sample.maternal_id]
  for pid in pids:
    if pid in L.sample_table:
      var p = L.sample_table[pid]
      if p.family_id == sample.family_id: continue
      # if p was already changed we have to use that id for the fam
      if sample.id in L.changed_samples:
        stderr.write_line &"[somalier] updating family_id for {p.id} to {sample.family_id}"
        p.family_id = sample.family_id
        L.changed_samples.incl(p.id)
        L.sample_table[p.id] = p
      elif p.id in L.changed_samples:
        sample.family_id = p.family_id
        stderr.write_line &"[somalier] updating family_id for {sample.id} to {sample.family_id}"
        L.sample_table[sample.id] = sample
        L.changed_samples.incl(sample.id)

      if p.family_id != sample.family_id:
        if L.sample_table[p.id].i < L.sample_table[sample.id].i:
          sample.family_id = p.family_id
          stderr.write_line &"[somalier] updating family_id for {sample.id} to {sample.family_id}"
          L.changed_samples.incl(sample.id)
          L.changed_samples.incl(sample.id)
        else:
          p.family_id = sample.family_id
          L.sample_table[p.id] = p
          stderr.write_line &"[somalier] updating family_id for {p.id} to {sample.family_id}"
          L.changed_samples.incl(p.id)


proc write_sample(fh:File, stats: seq[Stat4], gt_counts: array[5, seq[uint16]], i:int, L:SampleLooker) =
    let sample_name = L.sample_names[i]
    var sample = L.sample_table.getOrDefault(sample_name, Sample(id: sample_name, family_id: sample_name, sex: -9, phenotype:"-9", maternal_id:"-9", paternal_id:"-9"))
    fh.write(&"{sample.family_id}\t{sample.id}\t{sample.paternal_id}\t{sample.maternal_id}\t{sample.sex}\t{sample.phenotype}\t")
    fh.write(&"{L.sample_sex.getOrDefault(sample.id, \"-9\")}\t")
    fh.write(&"{stats[i].gtdp.mean:.1f}\t{stats[i].gtdp.standard_deviation():.1f}\t")
    fh.write(&"{stats[i].dp.mean:.1f}\t{stats[i].dp.standard_deviation():.1f}\t")
    fh.write(&"{stats[i].ab.mean:.2f}\t{stats[i].ab.standard_deviation():.2f}\t{gt_counts[0][i]}\t{gt_counts[1][i]}\t{gt_counts[2][i]}\t{gt_counts[3][i]}\t")
    fh.write(&"{gt_counts[4][i].float / (gt_counts[0][i] + gt_counts[1][i] + gt_counts[2][i] + gt_counts[3][i] + gt_counts[4][i]).float:.3f}\t")
    fh.write(&"{stats[i].x_dp.mean:.2f}\t{stats[i].x_dp.n}\t{stats[i].x_hom_ref}\t{stats[i].x_het}\t{stats[i].x_hom_alt}\t")
    fh.write(&"{stats[i].y_dp.mean:.2f}\t{stats[i].y_dp.n}\n")

proc look(final:relation_matrices, samples: var seq[Sample], stats: seq[Stat4], pairs: TableRef[string, seq[string]], sib_pairs: TableRef[string, seq[string]], relGt0p2: TableRef[string, seq[string]]): SampleLooker =
  result.sample_names = final.samples
  result.sample_table = newTable[string, Sample]()
  var tmp_sample_i = newTable[string, int]()
  for i, s in result.sample_names:
    tmp_sample_i[s] = i

  for s in samples.mitems:
    doAssert s.id notin result.sample_table, "error repeated sample id:" & s.id
    s.i = tmp_sample_i.getOrDefault(s.id, -1)
    result.sample_table[s.id] = s


  var byFam = newTable[string, seq[Sample]]()
  for i, s in samples.mpairs:
    byFam.mgetOrPut(s.family_id, @[]).add(s)

  for s in result.sample_names:
    if s notin result.sample_table:
      result.sample_table[s] = Sample(family_id: s, id: s, sex: -9, phenotype:"-9", maternal_id:"-9", paternal_id:"-9", i: result.sample_table.len)
      byFam.mgetOrPut(s, @[]).add(result.sample_table[s])

  # make sure sibs can join a family
  for a, bs in sib_pairs:
    var akid = result.sample_table[a]
    for b in bs:
      var bkid = result.sample_table[b]
      if akid.family_id != bkid.family_id:
        var bsamples: seq[Sample]
        doAssert byFam.take(bkid.family_id, bsamples)
        for bsample in bsamples:
          bsample.family_id = akid.family_id
          byFam[akid.family_id].add(bsample)
          result.sample_table[bsample.id].family_id = akid.family_id
  # join families on parent-child pairs
  for a, bs in pairs:
    var asample = result.sample_table[a]
    for b in bs:
      var bsamp = result.sample_table[b]
      if asample.family_id != bsamp.family_id:
        var bsamples: seq[Sample]
        doAssert byFam.take(bsamp.family_id, bsamples)
        for bsample in bsamples:
          bsample.family_id = asample.family_id
          byFam[asample.family_id].add(bsample)
          result.sample_table[bsample.id].family_id = asample.family_id

  # now reset family table with potentially updated ids
  byFam = newTable[string, seq[Sample]]()
  for s in result.sample_table.values:
    byFam.mgetOrPut(s.family_id, @[]).add(s)

  # join families with rel > 0.2 between any pair of samples
  for a, bs in relGt0p2:
    var asample = result.sample_table[a]
    for b in bs:
      var bsample = result.sample_table[b]
      if asample.family_id != bsample.family_id:
        stderr.write_line &"[somalier] joining families {asample.family_id} and {bsample.family_id} because of relatedness > 0.2"
        var bsamples: seq[Sample]
        doAssert byFam.take(bsample.family_id, bsamples)
        for bsample in bsamples:
          bsample.family_id = asample.family_id
          byFam[asample.family_id].add(bsample)
          result.sample_table[bsample.id].family_id = asample.family_id

  result.sample_sex = samples.to_sex_lookup
  result.has_y = stats.samples_have_y_depth
  result.pairs = pairs
  result.sib_pairs = sib_pairs


proc write_ped(fh:File, final: var relation_matrices, stats: seq[Stat4], gt_counts: array[5, seq[uint16]], L:var SampleLooker, infer:bool) =
  #var L = final.look(samples, stats, parent_child_pairs, sib_pairs)
  fh.write("#family_id\tsample_id\tpaternal_id\tmaternal_id\tsex\tphenotype\t")
  fh.write("original_pedigree_sex\tgt_depth_mean\tgt_depth_sd\tdepth_mean\tdepth_sd\tab_mean\tab_std\tn_hom_ref\tn_het\tn_hom_alt\tn_unknown\tp_middling_ab\t")
  fh.write("X_depth_mean\tX_n\tX_hom_ref\tX_het\tX_hom_alt\t")
  fh.write("Y_depth_mean\tY_n\n")

  if infer:
    final.remove_spurious_parent_ids(L, stats)
    for i, sample_name in L.sample_names:
      add_parents_and_check_sex(final, stats, gt_counts, i, L)
    add_siblings(final, stats, gt_counts, L)

    add_parent_to_sibs(final, stats, gt_counts, L)

    for i, sample_name in L.sample_names:
      update_family_ids(final, stats, gt_counts, i, L)
    final.remove_spurious_parent_ids(L, stats)

  for i, sample_name in L.sample_names:
    fh.write_sample(stats, gt_counts, i, L)
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

proc add_prefixed_samples(groups: var seq[pair], samples: seq[string], prefixes: seq[string]) =
  # update groups so that sample == ${prefix}sample
  #if len(prefixes) == 0: return
  let stripped = newTable[string, HashSet[string]]()

  for sample in samples:
    var s = sample
    for p in prefixes:
      if s.startsWith(p):
        s = s[p.len..s.high]
        break
    stripped.mgetOrPut(s, initHashSet[string](2)).incl(sample)

  for k, names in stripped:
    var names = names.toSeq
    if names.len < 2: continue
    for i in 0..<names.high:
      let A = names[i]
      for j in 1..names.high:
        let B = names[j]
        if A < B:
          groups.add((A, B, 1.0))
        else:
          groups.add((B, A, 1.0))

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
    option("--sample-prefix", multiple=true, help="optional sample prefixes that can be removed to find identical samples. e.g. batch1-sampleA batch2-sampleA")
    option("-p", "--ped", help="optional path to a ped/fam file indicating the expected relationships among samples.")
    option("-d", "--min-depth", default="7", help="only genotype sites with at least this depth.")
    option("--min-ab", default="0.3", help="hets sites must be between min-ab and 1 - min_ab. set this to 0.2 for RNA-Seq data")
    flag("-u", "--unknown", help="set unknown genotypes to hom-ref. it is often preferable to use this with VCF samples that were not jointly called")
    flag("-i", "--infer", help="infer relationships (https://github.com/brentp/somalier/wiki/pedigree-inference)")
    option("-o", "--output-prefix", help="output prefix for results.", default="somalier")
    arg("extracted", nargs= -1, help="$sample.somalier files for each sample. the first 10 are tested as a glob patterns")


  var opts = p.parse(argv)
  if opts.help:
    quit 0
  # first given 10 "files" could be a glob.
  opts.extracted.update_with_glob

  stderr.write_line &"[somalier] starting read of {opts.extracted.len} samples"
  if opts.extracted.len == 0 or (opts.extracted.len == 1 and not existsFile(opts.extracted[0])):
    echo p.help
    quit "[somalier] specify at least 1 extracted somalier file"
  var
    groups: seq[pair]
    samples: seq[Sample]
    min_depth = parseInt(opts.min_depth)
    min_ab = parseFloat(opts.min_ab)
    unk2hr = opts.unknown

  if not opts.output_prefix.endswith(".") or opts.output_prefix.endswith("/"):
    opts.output_prefix &= '.'

  var t0 = cpuTime()
  var final = read_extracted(opts.extracted, min_ab, min_depth, unk2hr)
  var n_samples = final.samples.len
  stderr.write_line &"[somalier] time to read files and get per-sample stats for {n_samples} samples: {cpuTime() - t0:.2f}"
  t0 = cpuTime()

  if opts.ped != "":
    samples = parse_ped(opts.ped)

  groups.add_ped_samples(samples, final.samples)

  groups.add_prefixed_samples(final.samples, opts.sample_prefix)
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
  fh_html.write(tmpls[0].replace("<SAMPLE_JSON>", toj(final.samples, final.stats, final.gt_counts, samples.to_sex_lookup)))

  var rels: seq[relations]

  var proportion_sampled = 200_000'f64 / float64(final.samples.len * final.samples.len)
  if proportion_sampled < 1:
    stderr.write_line &"[somalier] html and text output will have unrelated sample-pairs subset to {100 * proportion_sampled:.2f}% of points"

  fh_tsv.write_line '#', header.replace("$", "")
  var npairs:int
  var parent_child_pair = newTable[string, seq[string]]()
  var sib_pairs = newTable[string, seq[string]]()
  # we merge familys where the relatedness is > 0.2
  var relGt0p2 = newTable[string, seq[string]]()
  var nrels:int
  sort(groups, cmp_pair)
  for T in final.relatedness(grouped):
    let rel = T.r
    npairs.inc
    var idx = groups.binarySearch((rel.sample_a, rel.sample_b, -1.0), cmp_pair)
    if idx == -1:
      idx = groups.binarySearch((rel.sample_b, rel.sample_a, -1.0), cmp_pair)
    let expected_relatedness = if idx == -1: -1'f else: groups[idx].rel
    let rr = rel.rel

    # PARENT_CHILD
    if rr > 0.4 and rr < 0.6 and rel.ibs0.float / rel.ibs2.float < 0.005:
      parent_child_pair.mgetOrPut(rel.sample_a, @[]).add(rel.sample_b)
      parent_child_pair.mgetOrPut(rel.sample_b, @[]).add(rel.sample_a)

    elif rr > 0.38 and rr < 0.62 and rel.ibs0.float / rel.ibs2.float > 0.015 and rel.ibs0.float / rel.ibs2.float < 0.052:
      sib_pairs.mgetOrPut(rel.sample_a, @[]).add(rel.sample_b)
      sib_pairs.mgetOrPut(rel.sample_b, @[]).add(rel.sample_a)

    elif rr > 0.96 and rel.ibs0.float / rel.ibs2.float < 0.005:
      # this should be rare so we do a linear search
      var sibs = false;
      for pair in groups:
        if (pair.a == rel.sample_a and pair.b == rel.sample_b) or (pair.b == rel.sample_a and pair.a == rel.sample_b):
          if pair.rel > 0.4:
            sibs = true
          break
      if not sibs:
        if opts.infer:
          stderr.write_line &"[somalier] apparent identical twins or sample duplicate found with {rel.sample_a} and {rel.sample_b} NOT assuming siblings"
      else:
        if opts.infer:
          stderr.write_line &"[somalier] apparent identical twins or sample duplicate found with {rel.sample_a} and {rel.sample_b} assuming siblings as these were specified as such in the pedigree file"
          sib_pairs.mgetOrPut(rel.sample_a, @[]).add(rel.sample_b)
          sib_pairs.mgetOrPut(rel.sample_b, @[]).add(rel.sample_a)
    elif rr > 0.2 and final.gt_counts.high_quality(T.i) and final.gt_counts.high_quality(T.j):
      relGt0p2.mgetOrPut(rel.sample_a, @[]).add(rel.sample_b)
      relGt0p2.mgetOrPut(rel.sample_b, @[]).add(rel.sample_a)

    let ra = rand(1'f32)
    let interesting = expected_relatedness != -1 or rr > 0.05
    if (ra > proportion_sampled) and not interesting:
      continue
    rels.add(rel, max(0, expected_relatedness))
    nrels.inc

    fh_tsv.write_line rel.tsv(expected_relatedness)
  stderr.write_line &"[somalier] time to calculate all vs all relatedness for all {npairs} combinations: {cpuTime() - t0:.2f}"

  fh_html.write(%* rels)
  fh_html.write(tmpls[1])
  fh_html.close()
  stderr.write_line(&"[somalier] wrote interactive HTML output for {nrels} pairs to: ",  opts.output_prefix & "html")

  var L = final.look(samples, final.stats, parent_child_pair, sib_pairs, relGt0p2)
  fh_samples.write_ped(final, final.stats, final.gt_counts, L, opts.infer)

  fh_tsv.close()
  grouped.write(opts.output_prefix)

  stderr.write_line("[somalier] wrote groups to: ",  opts.output_prefix & "groups.tsv (look at this for cancer samples)")
  stderr.write_line("[somalier] wrote samples to: ",  opts.output_prefix & "samples.tsv")
  stderr.write_line("[somalier] wrote pair-wise relatedness metrics to: ",  opts.output_prefix & "pairs.tsv")
