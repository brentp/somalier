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

suite "CHARR estimator":
  test "hom-ref sites use stored alt AF in denominator":
    let stats = estimate_charr(@[ac(199, 1, 0)], @[0.25'f32], 7, 0.15)
    check abs(stats.contamination - 0.02) < 1e-8
    check stats.n_sites_usable == 1
    check stats.n_hom_ref_usable == 1
    check stats.n_hom_alt_usable == 0

  test "hom-alt sites use one minus stored alt AF in denominator":
    let stats = estimate_charr(@[ac(1, 199, 0)], @[0.2'f32], 7, 0.15)
    check abs(stats.contamination - 0.00625) < 1e-8
    check stats.n_sites_usable == 1
    check stats.n_hom_ref_usable == 0
    check stats.n_hom_alt_usable == 1

  test "looser CHARR cutoff keeps mildly contaminated homozygous sites":
    let strict = estimate_charr(@[ac(90, 10, 0)], @[0.4'f32], 7, 0.05)
    let loose = estimate_charr(@[ac(90, 10, 0)], @[0.4'f32], 7, 0.15)
    check strict.contamination.classify == fcNan
    check loose.n_sites_usable == 1
    check abs(loose.contamination - 0.25) < 1e-8

  test "het, low-depth, and high-other sites are excluded":
    let stats = estimate_charr(
      @[ac(199, 1, 0), ac(50, 50, 0), ac(6, 0, 0), ac(100, 0, 5)],
      @[0.25'f32, 0.25'f32, 0.25'f32, 0.25'f32],
      7,
      0.15
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
      0.15
    )
    check stats.contamination.classify == fcNan
    check stats.n_sites_usable == 0

suite "CHARR AF loading":
  test "flipped sites return AF aligned to stored B allele":
    let sites_path = "tests/_contam_flip_sites.vcf"
    write_sites(sites_path, @["1\t101\t.\tT\tA\t.\tPASS\tAF=0.2"])
    defer:
      if fileExists(sites_path):
        removeFile(sites_path)

    let site_data = readSitesWithAF(sites_path)
    check site_data.sites.len == 1
    check abs(site_data.pop_afs[0] - 0.8'f32) < 1e-6

suite "CHARR rows":
  test "sample-wise rows are emitted once per sketch":
    let
      workdir = "tests/_contam_rows"
      sites_path = workdir / "sites.vcf"
      out_path = workdir / "out.tsv"
      s1_path = workdir / "sample1.somalier"
      s2_path = workdir / "sample2.somalier"
    createDir(workdir)
    defer:
      if fileExists(out_path):
        removeFile(out_path)
      if fileExists(s1_path):
        removeFile(s1_path)
      if fileExists(s2_path):
        removeFile(s2_path)
      if fileExists(sites_path):
        removeFile(sites_path)
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

    let rows = charr_rows(@[s1_path, s2_path], sites_path, 7, 0.15)
    check rows.len == 2
    check rows[0].sample_name == "sample1"
    check rows[1].sample_name == "sample2"
    check abs(rows[0].stats.contamination - 0.01625) < 1e-8
    check abs(rows[1].stats.contamination - 0.0) < 1e-8

    write_charr_rows(rows, out_path)
    let lines = readFile(out_path).strip().splitLines()
    check lines.len == 3
    check lines[0] == "#sample_name\tn_sites_usable\tcontamination_charr"
    check lines[1].startsWith("sample1\t2\t")
    check lines[2].startsWith("sample2\t1\t")
