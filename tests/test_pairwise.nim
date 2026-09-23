import unittest

import somalierpkg/bitset
import somalierpkg/pairwise

type sample_data = object
  autosomal, x: genotypes
  gt_counts: array[5, uint16]

proc sample(calls: openArray[int8]): sample_data =
  result.autosomal.hom_ref = create_bitset(calls.len)
  result.autosomal.het = create_bitset(calls.len)
  result.autosomal.hom_alt = create_bitset(calls.len)
  result.x.hom_ref = create_bitset(1)
  result.x.het = create_bitset(1)
  result.x.hom_alt = create_bitset(1)
  for site, genotype in calls:
    if genotype >= 0:
      result.gt_counts[genotype].inc
      case genotype
      of 0: result.autosomal.hom_ref.set(site)
      of 1: result.autosomal.het.set(site)
      of 2: result.autosomal.hom_alt.set(site)
      else: discard
    else:
      result.gt_counts[3].inc

proc pair_view(value: sample_data): pair_sample =
  result.autosomal = value.autosomal.to_genotype_view
  result.x = value.x.to_genotype_view
  result.gt_counts = value.gt_counts

suite "matrix-free pair scoring":
  test "reports the established IBS and relatedness fields":
    let
      first = sample([0'i8, 1, 2, -1])
      second = sample([2'i8, 1, 2, 0])
      score = score_pair(first.pair_view, second.pair_view)
    check score.ibs0 == 1
    check score.ibs2 == 2
    check score.shared_hets == 1
    check score.shared_hom_alts == 1
    check score.het_ab == 2
    check score.n == 3
    check score.relatedness == -1.0
    check score.concordance == 0.0

  test "preserves the historical uint16 het-ab clamp":
    const site_count = uint16.high.int
    var calls = newSeq[int8](site_count)
    for genotype in calls.mitems:
      genotype = 1
    let value = sample(calls)
    let score = score_pair(value.pair_view, value.pair_view)
    check score.shared_hets == uint16.high
    check score.het_ab == uint16.high
    check score.relatedness == 2.0
