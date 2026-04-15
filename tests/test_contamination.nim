include somalierpkg/contamination

import unittest

proc ac(nref, nalt, nother: uint32): allele_count =
  allele_count(nref: nref, nalt: nalt, nother: nother)

proc write_sites(path: string, rows: openArray[string]) =
  var fh: File
  doAssert open(fh, path, fmWrite)
  fh.write_line("##fileformat=VCFv4.2")
  fh.write_line("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
  for row in rows:
    fh.write_line(row)
  fh.close()

proc counts_from_frac(p: float64, dp: int, nother: uint32 = 0'u32): allele_count =
  let nalt = uint32(int(round(p * dp.float64)))
  allele_count(nref: uint32(dp) - nalt, nalt: nalt, nother: nother)

suite "CHARR estimator":
  test "hom-ref sites use stored alt AF in denominator":
    let stats = estimate_charr(@[ac(199, 1, 0)], @[0.25'f32], 7, 0.10, 0.001)
    check abs(stats.contamination - 0.02) < 1e-8
    check stats.n_sites_usable == 1
    check stats.n_hom_ref_usable == 1
    check stats.n_hom_alt_usable == 0

  test "hom-alt sites use one minus stored alt AF in denominator":
    let stats = estimate_charr(@[ac(1, 199, 0)], @[0.2'f32], 7, 0.10, 0.001)
    check abs(stats.contamination - 0.00625) < 1e-8
    check stats.n_sites_usable == 1
    check stats.n_hom_ref_usable == 0
    check stats.n_hom_alt_usable == 1

  test "depth-aware CHARR rule keeps moderate imbalance but rejects extreme high-depth imbalance":
    let moderate = estimate_charr(@[ac(90, 10, 0)], @[0.4'f32], 7, 0.10, 0.001)
    let extreme = estimate_charr(@[ac(60, 40, 0)], @[0.4'f32], 7, 0.10, 0.001)
    check moderate.n_sites_usable == 1
    check abs(moderate.contamination - 0.25) < 1e-8
    check extreme.contamination.classify == fcNan

  test "het, low-depth, and high-other sites are excluded":
    let stats = estimate_charr(
      @[ac(199, 1, 0), ac(50, 50, 0), ac(6, 0, 0), ac(100, 0, 5)],
      @[0.25'f32, 0.25'f32, 0.25'f32, 0.25'f32],
      7,
      0.10,
      0.001
    )
    check abs(stats.contamination - 0.02) < 1e-8
    check stats.n_sites_usable == 1
    check stats.n_hom_ref_usable == 1
    check stats.n_hom_alt_usable == 0

  test "zero usable sites returns NaN":
    let stats = estimate_charr(
      @[ac(50, 50, 0), ac(6, 0, 0), ac(100, 0, 5)],
      @[0.25'f32, 0.25'f32, 0.25'f32],
      7,
      0.10,
      0.001
    )
    check stats.contamination.classify == fcNan
    check stats.n_sites_usable == 0

suite "Exact pairwise helpers":
  test "clean hom-ref plus hom-alt contaminant respects alpha and error":
    check abs(observed_alt_prob(0'i8, 2'i8, 0.2, 0.0) - 0.2) < 1e-12
    check abs(observed_alt_prob(0'i8, 2'i8, 0.2, 0.002) - 0.2012) < 1e-12

  test "clean hom-alt plus hom-ref contaminant respects alpha and error":
    check abs(observed_alt_prob(2'i8, 0'i8, 0.2, 0.0) - 0.8) < 1e-12
    check abs(observed_alt_prob(2'i8, 0'i8, 0.2, 0.002) - 0.7988) < 1e-12

  test "site_log_likelihood clamps probabilities away from zero and one":
    let hom_ref = site_log_likelihood(100, 0, 0.0)
    let hom_alt = site_log_likelihood(0, 100, 1.0)
    check hom_ref.classify == fcNormal
    check hom_alt.classify == fcNormal

suite "Exact pairwise estimator":
  test "same-sample-style anchor no longer lands at the boundary":
    let
      pop_afs = @[0.2'f32, 0.8'f32, 0.3'f32, 0.7'f32]
      sample_a = @[ac(200, 0, 0), ac(0, 200, 0), ac(200, 0, 0), ac(0, 200, 0)]
      a_from_a = estimate_pair_contamination(sample_a, sample_a, pop_afs, 7)
    check a_from_a.n_sites_usable == 4
    check a_from_a.contamination < 0.02
    check a_from_a.contamination < 0.5

  test "synthetic contamination is directional with the correct anchor":
    let
      pop_afs = @[0.25'f32, 0.25'f32, 0.75'f32, 0.75'f32]
      anchor = @[
        counts_from_frac(0.0, 200),
        counts_from_frac(0.0, 200),
        counts_from_frac(1.0, 200),
        counts_from_frac(1.0, 200)
      ]
      receiver = @[
        counts_from_frac(0.2, 200),
        counts_from_frac(0.2, 200),
        counts_from_frac(0.8, 200),
        counts_from_frac(0.8, 200)
      ]
      receiver_from_anchor = estimate_pair_contamination(receiver, anchor, pop_afs, 7)
      anchor_from_receiver = estimate_pair_contamination(anchor, receiver, pop_afs, 7)
    check receiver_from_anchor.n_sites_usable == pop_afs.len
    check receiver_from_anchor.contamination > 0.10
    check anchor_from_receiver.contamination.classify == fcNan or
        anchor_from_receiver.contamination + 0.03 < receiver_from_anchor.contamination

  test "zero usable sites returns NaN":
    let stats = estimate_pair_contamination(
      @[ac(6, 0, 0)],
      @[ac(0, 6, 0)],
      @[0.2'f32],
      7
    )
    check stats.n_sites_usable == 0
    check stats.contamination.classify == fcNan

  test "anchor het, low-depth, and high-other sites are excluded":
    let stats = estimate_pair_contamination(
      @[ac(100, 0, 0), ac(100, 0, 0), ac(100, 0, 0), ac(100, 0, 0)],
      @[ac(0, 100, 0), ac(50, 50, 0), ac(6, 0, 0), ac(100, 0, 5)],
      @[0.2'f32, 0.2'f32, 0.2'f32, 0.2'f32],
      7
    )
    check stats.n_sites_usable == 1

suite "AF loading":
  test "flipped sites return AF aligned to stored B allele":
    let sites_path = "tests/_contam_flip_sites.vcf"
    write_sites(sites_path, @["1\t101\t.\tT\tA\t.\tPASS\tAF=0.2"])
    defer:
      if fileExists(sites_path):
        removeFile(sites_path)

    let site_data = readSitesWithAF(sites_path)
    check site_data.sites.len == 1
    check abs(site_data.pop_afs[0] - 0.8'f32) < 1e-6

suite "Contamination outputs":
  test "sample and pair rows are emitted from shared loaded sketches":
    let
      workdir = "tests/_contam_rows"
      sites_path = workdir / "sites.vcf"
      sample_prefix = workdir / "cont"
      sample_path = sample_output_path(sample_prefix)
      pair_path = pair_output_path(sample_prefix)
      s1_path = workdir / "sample1.somalier"
      s2_path = workdir / "sample2.somalier"
    createDir(workdir)
    defer:
      for path in [pair_path, sample_path, s1_path, s2_path, sites_path]:
        if fileExists(path):
          removeFile(path)
      if dirExists(workdir):
        removeDir(workdir)

    write_sites(sites_path, @[
      "1\t101\t.\tA\tC\t.\tPASS\tAF=0.25",
      "1\t201\t.\tG\tT\t.\tPASS\tAF=0.4"
    ])

    write_counts(counts(sample_name: "sample1",
      sites: @[ac(199, 1, 0), ac(199, 1, 0)],
      x_sites: @[],
      y_sites: @[]), "sample1", s1_path)
    write_counts(counts(sample_name: "sample2",
      sites: @[ac(100, 0, 0), ac(80, 20, 0)],
      x_sites: @[],
      y_sites: @[]), "sample2", s2_path)

    let loaded = load_contamination_inputs(@[s1_path, s2_path], sites_path)
    let sample_rows = charr_rows(loaded.sketches, loaded.pop_afs, 7, 0.10, 0.001)
    let pairs = pair_rows(loaded.sketches, loaded.log_priors, defaultExactErrorRate)
    check loaded.log_priors.len == 2
    check abs(loaded.log_priors[0][0] - ln(0.75 * 0.75)) < 1e-12
    check abs(loaded.log_priors[0][1] - ln(2.0 * 0.25 * 0.75)) < 1e-12
    check abs(loaded.log_priors[0][2] - ln(0.25 * 0.25)) < 1e-12
    check sample_rows.len == 2
    check sample_rows[0].sample_name == "sample1"
    check sample_rows[1].sample_name == "sample2"
    check abs(sample_rows[0].stats.contamination - 0.01625) < 1e-8
    check abs(sample_rows[1].stats.contamination - 0.25) < 1e-8

    check pairs.len == 2
    check pairs[0].sample_name == "sample1"
    check pairs[0].anchor_sample == "sample2"
    check pairs[1].sample_name == "sample2"
    check pairs[1].anchor_sample == "sample1"

    write_charr_rows(sample_rows, sample_path)
    write_pair_rows(pairs, pair_path)

    let sample_lines = readFile(sample_path).strip().splitLines()
    check sample_lines.len == 3
    check sample_lines[0] == "#sample_name\tn_sites_usable\tcontamination_charr"
    check sample_lines[1].startsWith("sample1\t2\t")
    check sample_lines[2].startsWith("sample2\t2\t")

    let pair_lines = readFile(pair_path).strip().splitLines()
    check pair_lines.len == 3
    check pair_lines[0] == "#sample_name\tanchor_sample\tn_sites_usable\tcontamination_mle"
    check pair_lines[1].startsWith("sample1\tsample2\t")
    check pair_lines[2].startsWith("sample2\tsample1\t")

  test "tumor-normal pair mode emits a single tumor-target normal-reference row":
    let
      workdir = "tests/_contam_pair_mode"
      sites_path = workdir / "sites.vcf"
      tumor_path = workdir / "tumor.somalier"
      normal_path = workdir / "normal.somalier"
    createDir(workdir)
    defer:
      for path in [tumor_path, normal_path, sites_path]:
        if fileExists(path):
          removeFile(path)
      if dirExists(workdir):
        removeDir(workdir)

    write_sites(sites_path, @[
      "1\t101\t.\tA\tC\t.\tPASS\tAF=0.25",
      "1\t201\t.\tG\tT\t.\tPASS\tAF=0.75"
    ])

    write_counts(counts(sample_name: "tumor",
      sites: @[ac(160, 40, 0), ac(40, 160, 0)],
      x_sites: @[],
      y_sites: @[]), "tumor", tumor_path)
    write_counts(counts(sample_name: "normal",
      sites: @[ac(200, 0, 0), ac(0, 200, 0)],
      x_sites: @[],
      y_sites: @[]), "normal", normal_path)

    let loaded = load_contamination_inputs(@[tumor_path, normal_path], sites_path)
    let row = tumor_normal_pair_row(loaded, defaultExactErrorRate)
    check row.sample_name == "tumor"
    check row.anchor_sample == "normal"
    check row.stats.n_sites_usable == 2
    check pair_row_line(row).startsWith("tumor\tnormal\t2\t")
    check row.stats.contamination > 0.10

  test "tumor-normal pair flag uses exactly two extracted inputs in order":
    let
      workdir = "tests/_contam_pair_args"
      tumor_path = workdir / "tumor.somalier"
      normal_path = workdir / "normal.somalier"
    createDir(workdir)
    defer:
      for path in [tumor_path, normal_path]:
        if fileExists(path):
          removeFile(path)
      if dirExists(workdir):
        removeDir(workdir)

    writeFile(tumor_path, "")
    writeFile(normal_path, "")
    let cfg = parse_contamination_args(@[
      "contamination",
      "-s", "tests/test_sites.vcf",
      "-p",
      tumor_path,
      normal_path
    ])
    check cfg.tumor_normal_pair_mode
    check cfg.extracted == @[tumor_path, normal_path]

  test "output prefix strips tsv suffix and falls back to default":
    check normalize_output_prefix("cont.tsv") == "cont"
    check sample_output_path(normalize_output_prefix("cont.tsv")) == "cont.samples.tsv"
    check pair_output_path(normalize_output_prefix("cont.tsv")) == "cont.pairs.tsv"
    check normalize_output_prefix("") == defaultOutputPrefix
