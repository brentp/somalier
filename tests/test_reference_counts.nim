include somalierpkg/relate

import unittest

proc sample_json(r: relation_matrices): JsonNode =
  parseJson(toj(r.samples, r.stats, r.gt_counts, @[CharrStats()],
      newTable[string, string]()))[0]

proc reported_counts(r: relation_matrices, prefix: string = "n_"): seq[int] =
  let sample = r.sample_json
  for genotype in ["hom_ref", "het", "hom_alt"]:
    result.add(sample[prefix & genotype].getInt)

suite "reference-oriented sample counts":
  let workdir = getTempDir() / ("somalier-reference-counts-" & $getCurrentProcessId())
  let sites_path = workdir / "sites.vcf"
  let sketch_path = workdir / "sample.somalier"

  setup:
    createDir(workdir)

  teardown:
    removeDir(workdir)

  test "existing sketches report 20 hom-ref while retaining A/B genotypes":
    var rows = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    var sample = counts(sample_name: "test_sample")
    # Reverse the sites file to exercise the same sorting used by extract.
    for i in countdown(19, 0):
      let alleles = if i mod 2 == 0: "G\tC" else: "A\tG"
      rows.add(&"1\t{i + 1}\t.\t{alleles}\t.\tPASS\tAF=0.2\n")
    for i in 0..<20:
      sample.sites.add(if i mod 2 == 0: allele_count(nalt: 80)
                      else: allele_count(nref: 80))
    sample.sites.add(allele_count())
    rows.add("1\t21\t.\tT\tC\t.\tPASS\tAF=0.2\n")
    writeFile(sites_path, rows)
    sample.write_counts(sample.sample_name, sketch_path)
    let original = readFile(sketch_path)

    let legacy = read_extracted(@[sketch_path], 0.3, 7, false)
    let oriented = read_extracted(@[sketch_path], 0.3, 7, false, sites_path)
    check legacy.reported_counts == @[10, 0, 10]
    check oriented.reported_counts == @[20, 0, 0]
    check oriented.gt_counts[3][0] == 1
    check oriented.gt_counts == legacy.gt_counts
    check oriented.genotypes == legacy.genotypes
    check oriented.allele_counts == legacy.allele_counts
    check readFile(sketch_path) == original
    check original[0].uint8 == 2'u8

    let encoded = oriented.sample_json
    check encoded["n_hom_ref"].getInt == 20
    check encoded["n_hom_alt"].getInt == 0
    check encoded["n_unknown"].getInt == 1

    let output_path = workdir / "samples.tsv"
    let output = open(output_path, fmWrite)
    let looker = SampleLooker(sample_names: oriented.samples,
        sample_table: newTable[string, Sample](), sample_sex: newTable[string, string]())
    output.write_sample(oriented.stats, oriented.gt_counts, @[CharrStats()], 0, looker)
    output.close()
    let fields = readFile(output_path).strip.split('\t')
    check fields[13..16] == @["20", "0", "0", "1"]

    let imputed = read_extracted(@[sketch_path], 0.3, 7, true, sites_path)
    check imputed.reported_counts == @[21, 0, 0]
    check imputed.gt_counts[3][0] == 0
    let legacy_unknown = read_extracted(@[sketch_path], 0.3, 7, true)
    check imputed.gt_counts == legacy_unknown.gt_counts
    check imputed.genotypes == legacy_unknown.genotypes

  test "hom-alt, het, unknown and X counts use the correct orientation":
    writeFile(sites_path, "1\t1\t.\tT\tA\n1\t2\t.\tG\tA\n" &
        "1\t3\t.\tC\tA\n1\t4\t.\tT\tA\n" &
        "chrX\t1\t.\tT\tC\nchrX\t2\t.\tG\tA\n" &
        "chrY\t1\t.\tT\tA\n")
    counts(sites: @[allele_count(nref: 80), allele_count(nref: 40, nalt: 40),
                   allele_count(nref: 1), allele_count(nref: 80, nalt: 20)],
           x_sites: @[allele_count(nalt: 80), allele_count(nref: 40, nalt: 40)],
           y_sites: @[allele_count(nalt: 80)]).write_counts("sample", sketch_path)
    let oriented = read_extracted(@[sketch_path], 0.3, 7, false, sites_path)
    check oriented.reported_counts == @[0, 1, 1]
    check oriented.gt_counts[3][0] == 2
    check oriented.reported_counts("x_") == @[1, 1, 0]
    let imputed = read_extracted(@[sketch_path], 0.3, 7, true, sites_path)
    check imputed.reported_counts == @[2, 1, 1]

  test "rejects mismatched autosomal, X and Y site counts":
    counts(sites: @[allele_count(nref: 80)],
           x_sites: @[allele_count(nref: 80)],
           y_sites: @[allele_count(nref: 80)]).write_counts("sample", sketch_path)
    let rows = @["1\t1\t.\tA\tG\n", "X\t1\t.\tA\tG\n", "Y\t1\t.\tA\tG\n"]
    for omitted in 0..2:
      var sites = ""
      for i, row in rows:
        if i != omitted: sites.add(row)
      writeFile(sites_path, sites)
      expect ValueError:
        discard read_extracted(@[sketch_path], 0.3, 7, false, sites_path)
