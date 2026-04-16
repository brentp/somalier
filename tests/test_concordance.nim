include somalierpkg/relate

import unittest

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

suite "inferred concordance":
  test "uses the smaller callable homozygous set as the denominator":
    var rm = relation_matrices(
      genotypes: @[
        build_genotypes([0], [2], [1], 4),
        build_genotypes([0, 2], [], [1, 3], 4)
      ],
      gt_counts: [
        @[1'u16, 0'u16],
        @[1'u16, 0'u16],
        @[1'u16, 2'u16],
        @[0'u16, 0'u16],
        @[0'u16, 0'u16],
      ]
    )
    check abs(rm.inferred_hom_concordance(0, 1) - 1.0'f32) < 0.0001

  test "ties still compute concordance from callable homozygous markers":
    var rm = relation_matrices(
      genotypes: @[
        build_genotypes([0], [2], [1], 4),
        build_genotypes([0, 2], [1], [3], 4)
      ],
      gt_counts: [
        @[1'u16, 1'u16],
        @[1'u16, 1'u16],
        @[1'u16, 1'u16],
        @[0'u16, 0'u16],
        @[0'u16, 0'u16],
      ]
    )

    check abs(rm.inferred_hom_concordance(0, 1) - 0.5'f32) < 0.0001

suite "adjusted concordance":
  test "clean pairs keep perfect concordance":
    var rm = relation_matrices(
      genotypes: @[
        build_genotypes([0, 1, 2, 3], [], [], 4),
        build_genotypes([0, 1, 2, 3], [], [], 4)
      ],
      gt_counts: [
        @[4'u16, 4'u16],
        @[0'u16, 0'u16],
        @[0'u16, 0'u16],
        @[0'u16, 0'u16],
        @[0'u16, 0'u16],
      ]
    )
    let rel = relation(shared_hom_alts: 2, ibs0: 0, hom_alts_a: 2, hom_alts_b: 2)

    check abs(rm.adjusted_concordance(rel, 0, 1) - 1.0'f32) < 0.0001

  test "middling allele balance increases the low hom-alt penalty":
    var clean = relation_matrices(
      genotypes: @[
        build_genotypes([0, 1, 2, 3], [], [], 4),
        build_genotypes([0, 1, 2, 3], [], [], 4)
      ],
      gt_counts: [
        @[4'u16, 4'u16],
        @[0'u16, 0'u16],
        @[0'u16, 0'u16],
        @[0'u16, 0'u16],
        @[0'u16, 0'u16],
      ]
    )
    var noisy = clean
    noisy.gt_counts[4] = @[1'u16, 1'u16]

    let rel = relation(shared_hom_alts: 4, ibs0: 0, hom_alts_a: 5, hom_alts_b: 5)

    let clean_score = clean.adjusted_concordance(rel, 0, 1)
    let noisy_score = noisy.adjusted_concordance(rel, 0, 1)

    check abs(raw_hom_alt_concordance(rel) - 0.8'f32) < 0.0001
    check clean_score > noisy_score
    check clean_score > 0.99'f32
    check noisy_score > 0.8'f32
