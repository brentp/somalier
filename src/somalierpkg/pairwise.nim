## Matrix-free exact Somalier pair scoring.
##
## Keep this module independent of candidate generation. Exhaustive and sparse
## callers must use the same scorer so retrieval can change which pairs are
## examined without changing any pair statistic.

import bitops
import math
import ./bitset

type
  ## Non-owning view over fixed genotype words. The underlying arena or bitset
  ## must remain allocated and stationary for the lifetime of this view.
  bitset_view* = object
    words*: ptr UncheckedArray[uint64]
    len*: int

  genotype_view* = object
    hom_ref*: bitset_view
    het*: bitset_view
    hom_alt*: bitset_view

  pair_sample* = object
    autosomal*: genotype_view
    x*: genotype_view
    gt_counts*: array[5, uint16]

  pair_score* = object
    concordance*: float32
    hets_a*: uint16
    hets_b*: uint16
    hom_alts_a*: uint16
    hom_alts_b*: uint16
    shared_hom_alts*: uint16
    shared_hets*: uint16
    het_ab*: uint16
    ibs0*: uint16
    ibs2*: uint16
    x_ibs0*: uint16
    x_ibs2*: uint16
    n*: uint16

proc to_bitset_view*(value: bitset): bitset_view {.inline.} =
  result.len = value.len
  if value.len > 0:
    result.words = cast[ptr UncheckedArray[uint64]](unsafeAddr value[0])

proc to_genotype_view*(value: genotypes): genotype_view {.inline.} =
  result.hom_ref = value.hom_ref.to_bitset_view
  result.het = value.het.to_bitset_view
  result.hom_alt = value.hom_alt.to_bitset_view

proc check_compatible(first, second: genotype_view) {.inline.} =
  doAssert first.hom_ref.len == first.het.len and
    first.hom_ref.len == first.hom_alt.len
  doAssert second.hom_ref.len == second.het.len and
    second.hom_ref.len == second.hom_alt.len
  doAssert first.hom_ref.len == second.hom_ref.len

proc relatedness*(score: pair_score): float64 {.inline.} =
  2 * (score.shared_hets.float64 - 2 * score.ibs0.float64) /
    max(1, score.het_ab.float64)

proc clamp_01(value: float32): float32 {.inline.} =
  max(0'f32, min(1'f32, value))

proc stretch_concordance(value: float32): float32 {.inline.} =
  const anchor = 0.4'f32
  if value <= anchor:
    return value
  let scaled = (value - anchor) / (1'f32 - anchor)
  anchor + (1'f32 - anchor) * cbrt(scaled)

proc p_middling_ab(sample: pair_sample): float32 {.inline.} =
  let total = max(1'u16, sample.gt_counts[0] + sample.gt_counts[1] +
      sample.gt_counts[2] + sample.gt_counts[3] + sample.gt_counts[4]).float32
  sample.gt_counts[4].float32 / total

proc inferred_hom_concordance*(first, second: pair_sample): float32 {.inline.} =
  check_compatible(first.autosomal, second.autosomal)
  var first_ref_sites = 0
  var second_ref_sites = 0
  var matches = 0
  for idx in 0 ..< first.autosomal.hom_ref.len:
    let first_known = first.autosomal.hom_ref.words[idx] or
      first.autosomal.het.words[idx] or first.autosomal.hom_alt.words[idx]
    let second_known = second.autosomal.hom_ref.words[idx] or
      second.autosomal.het.words[idx] or second.autosomal.hom_alt.words[idx]
    let first_hom = first.autosomal.hom_ref.words[idx] or
      first.autosomal.hom_alt.words[idx]
    let second_hom = second.autosomal.hom_ref.words[idx] or
      second.autosomal.hom_alt.words[idx]
    first_ref_sites += countSetBits(first_hom and second_known).int
    second_ref_sites += countSetBits(second_hom and first_known).int
    matches += countSetBits(
      (first.autosomal.hom_ref.words[idx] and second.autosomal.hom_ref.words[idx]) or
      (first.autosomal.hom_alt.words[idx] and second.autosomal.hom_alt.words[idx])).int
  let denominator = max(1, min(first_ref_sites, second_ref_sites))
  matches.float32 / denominator.float32

proc raw_hom_alt_concordance*(score: pair_score): float32 {.inline.} =
  (score.shared_hom_alts.float32 - 2'f32 * score.ibs0.float32) /
    max(1'u16, min(score.hom_alts_a, score.hom_alts_b)).float32

proc adjusted_concordance*(first, second: pair_sample;
                          score: pair_score): float32 {.inline.} =
  let base = inferred_hom_concordance(first, second)
  let hom_alt = clamp_01(score.raw_hom_alt_concordance)
  let pm = (first.p_middling_ab + second.p_middling_ab) / 2'f32
  let low_hom_alt_penalty = max(0'f32, 0.70'f32 - hom_alt + 2'f32 * pm)
  stretch_concordance(clamp_01(base - low_hom_alt_penalty))

proc score_pair*(first, second: pair_sample): pair_score =
  check_compatible(first.autosomal, second.autosomal)
  check_compatible(first.x, second.x)
  result = pair_score(
    hets_a: first.gt_counts[1],
    hets_b: second.gt_counts[1],
    hom_alts_a: first.gt_counts[2],
    hom_alts_b: second.gt_counts[2],
  )
  var het_ab: int32
  for idx in 0 ..< first.autosomal.hom_ref.len:
    let
      first_hom_ref = first.autosomal.hom_ref.words[idx]
      first_het = first.autosomal.het.words[idx]
      first_hom_alt = first.autosomal.hom_alt.words[idx]
      second_hom_ref = second.autosomal.hom_ref.words[idx]
      second_het = second.autosomal.het.words[idx]
      second_hom_alt = second.autosomal.hom_alt.words[idx]
      jointly_called = (first_hom_ref or first_het or first_hom_alt) and
        (second_hom_ref or second_het or second_hom_alt)
      shared_hets = first_het and second_het
      shared_hom_alts = first_hom_alt and second_hom_alt
    result.ibs0 += countSetBits(
      (first_hom_ref and second_hom_alt) or
      (first_hom_alt and second_hom_ref)).uint16
    result.ibs2 += countSetBits(shared_hom_alts or shared_hets or
      (first_hom_ref and second_hom_ref)).uint16
    result.shared_hets += countSetBits(shared_hets).uint16
    result.shared_hom_alts += countSetBits(shared_hom_alts).uint16
    het_ab += (countSetBits(second_het and jointly_called) +
      countSetBits(first_het and jointly_called)).int32
    result.n += countSetBits(jointly_called).uint16
  result.het_ab = min(uint16.high.int32, het_ab).uint16
  for idx in 0 ..< first.x.hom_ref.len:
    result.x_ibs0 += countSetBits(
      (first.x.hom_ref.words[idx] and second.x.hom_alt.words[idx]) or
      (first.x.hom_alt.words[idx] and second.x.hom_ref.words[idx])).uint16
    result.x_ibs2 += countSetBits(
      (first.x.hom_alt.words[idx] and second.x.hom_alt.words[idx]) or
      (first.x.het.words[idx] and second.x.het.words[idx]) or
      (first.x.hom_ref.words[idx] and second.x.hom_ref.words[idx])).uint16
  result.concordance = adjusted_concordance(first, second, result)
