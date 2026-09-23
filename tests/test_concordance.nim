import unittest
import somalierpkg/bitset
import somalierpkg/pairwise

type sample_data = object
  autosomal: genotypes
  gt_counts: array[5, uint16]

proc build_genotypes(hom_refs: openArray[int], hets: openArray[int],
    hom_alts: openArray[int], n_sites: int): genotypes =
  result.hom_ref = create_bitset(n_sites)
  result.het = create_bitset(n_sites)
  result.hom_alt = create_bitset(n_sites)
  for i in hom_refs:
    result.hom_ref.set(i)
  for i in hets:
    result.het.set(i)
  for i in hom_alts:
    result.hom_alt.set(i)

proc build_sample(hom_refs: openArray[int], hets: openArray[int],
    hom_alts: openArray[int], n_sites: int,
    gt_counts: array[5, uint16]): sample_data =
  result.autosomal = build_genotypes(hom_refs, hets, hom_alts, n_sites)
  result.gt_counts = gt_counts

proc pair_view(value: sample_data): pair_sample =
  result.autosomal = value.autosomal.to_genotype_view
  result.gt_counts = value.gt_counts

suite "inferred concordance":
  test "uses the smaller callable homozygous set as the denominator":
    let first = build_sample([0], [2], [1], 4,
      [1'u16, 1'u16, 1'u16, 0'u16, 0'u16])
    let second = build_sample([0, 2], [], [1, 3], 4,
      [0'u16, 0'u16, 2'u16, 0'u16, 0'u16])
    check abs(inferred_hom_concordance(first.pair_view,
      second.pair_view) - 1.0'f32) < 0.0001

  test "ties still compute concordance from callable homozygous markers":
    let first = build_sample([0], [2], [1], 4,
      [1'u16, 1'u16, 1'u16, 0'u16, 0'u16])
    let second = build_sample([0, 2], [1], [3], 4,
      [1'u16, 1'u16, 1'u16, 0'u16, 0'u16])
    check abs(inferred_hom_concordance(first.pair_view,
      second.pair_view) - 0.5'f32) < 0.0001

suite "adjusted concordance":
  test "clean pairs keep perfect concordance":
    let first = build_sample([0, 1, 2, 3], [], [], 4,
      [4'u16, 0'u16, 0'u16, 0'u16, 0'u16])
    let second = first
    let score = pair_score(shared_hom_alts: 2, ibs0: 0,
      hom_alts_a: 2, hom_alts_b: 2)
    check abs(adjusted_concordance(first.pair_view,
      second.pair_view, score) - 1.0'f32) < 0.0001

  test "middling allele balance increases the low hom-alt penalty":
    let clean = build_sample([0, 1, 2, 3], [], [], 4,
      [4'u16, 0'u16, 0'u16, 0'u16, 0'u16])
    var noisy = clean
    noisy.gt_counts[4] = 1
    let score = pair_score(shared_hom_alts: 4, ibs0: 0,
      hom_alts_a: 5, hom_alts_b: 5)

    let clean_score = adjusted_concordance(clean.pair_view,
      clean.pair_view, score)
    let noisy_score = adjusted_concordance(noisy.pair_view,
      noisy.pair_view, score)

    check abs(raw_hom_alt_concordance(score) - 0.8'f32) < 0.0001
    check clean_score > noisy_score
    check clean_score > 0.99'f32
    check noisy_score > 0.8'f32
