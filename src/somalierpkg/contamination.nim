import os
import streams
import strutils
import tables
import strformat
import math
import argparse
import ./common
import ./depthview
import ./relate
import ./pedrel

type ContaminationConfig = object
  sites: string
  min_depth: int
  charr_hom_cutoff: float
  output: string
  extracted: seq[string]

type CharrStats* = object
  contamination*: float64
  n_sites_usable*: int
  n_hom_ref_usable*: int
  n_hom_alt_usable*: int

type CharrRow = tuple[sample_name: string, stats: CharrStats]

type SketchMeta = object
  path: string
  sample_name: string
  n_sites: uint16
  nx_sites: uint16
  ny_sites: uint16

proc read_bytes(f: FileStream, p: pointer, n: int, path: string, field: string) =
  if n == 0:
    return
  if n != f.readData(p, n):
    quit &"[somalier] unexpected end of file while reading {field} from {path}"

proc read_sketch_header(path: string): SketchMeta =
  var f = newFileStream(path, fmRead)
  if f == nil:
    quit "[somalier] could not open sketch file: " & path
  defer:
    f.close()

  var version: uint8
  read_bytes(f, version.addr, sizeof(version), path, "format version")
  doAssert version == formatVersion, &"expected matching versions got {version}, expected {formatVersion}"

  var sample_len: uint8
  read_bytes(f, sample_len.addr, sizeof(sample_len), path, "sample name length")
  result.sample_name = newString(sample_len)
  if sample_len > 0'u8:
    read_bytes(f, result.sample_name[0].addr, sample_len.int, path, "sample name")

  read_bytes(f, result.n_sites.addr, sizeof(result.n_sites), path, "autosomal site count")
  read_bytes(f, result.nx_sites.addr, sizeof(result.nx_sites), path, "x site count")
  read_bytes(f, result.ny_sites.addr, sizeof(result.ny_sites), path, "y site count")
  result.path = path

proc read_sketch_headers(paths: seq[string], expected_autosomal_sites: int): seq[SketchMeta] =
  var seen = initTable[string, bool]()
  for path in paths:
    let meta = read_sketch_header(path)
    if meta.sample_name.len == 0:
      quit "[somalier] empty sample name found in sketch file: " & path
    if meta.sample_name in seen:
      quit "[somalier] duplicate sample name '" & meta.sample_name & "' found in extracted inputs"
    seen[meta.sample_name] = true
    # Version 1 relies on the user giving the same sites file used at extract
    # time, so fail early on any autosomal sketch-length mismatch.
    if meta.n_sites.int != expected_autosomal_sites:
      quit &"[somalier] sketch file '{path}' has {meta.n_sites} autosomal sites, expected {expected_autosomal_sites} from --sites"
    result.add(meta)

proc parse_contamination_args(argv: seq[string]): ContaminationConfig =
  var argv = argv
  if argv.len > 0 and argv[0] == "contamination":
    argv = argv[1..argv.high]
  if argv.len == 0:
    argv = @["-h"]

  var p = newParser("somalier contamination"):
    help("estimate contamination from extracted somalier files")
    option("-s", "--sites",
      help = "sites VCF with AF in INFO; this must match the site set used during extract")
    option("--min-depth", default = "7",
      help = "minimum depth required for a site to be usable")
    option("--charr-hom-cutoff", default = "0.15",
      help = "maximum minor-allele fraction still treated as homozygous for the CHARR summary estimate")
    option("-o", "--output", default = "-",
      help = "output TSV path (default: stdout)")
    arg("extracted", nargs = -1,
      help = "$sample.somalier files for each sample. the first 10 are tested as glob patterns")

  let opts = p.parse(argv)
  if opts.help:
    quit 0

  if opts.sites == "":
    echo p.help
    quit "[somalier] --sites file required"
  if not fileExists(opts.sites):
    quit "[somalier] unable to open sites file: " & opts.sites

  result.extracted = opts.extracted
  result.extracted.update_with_glob
  if result.extracted.len == 0 or
      (result.extracted.len == 1 and not fileExists(result.extracted[0])):
    echo p.help
    quit "[somalier] specify at least 1 extracted somalier file"

  try:
    result.min_depth = parseInt(opts.min_depth)
  except ValueError:
    quit "[somalier] --min-depth must be an integer"
  if result.min_depth < 1:
    quit "[somalier] --min-depth must be at least 1"
  try:
    result.charr_hom_cutoff = parseFloat(opts.charr_hom_cutoff)
  except ValueError:
    quit "[somalier] --charr-hom-cutoff must be a number"
  if result.charr_hom_cutoff <= 0 or result.charr_hom_cutoff >= 0.5:
    quit "[somalier] --charr-hom-cutoff must be greater than 0 and less than 0.5"

  result.sites = opts.sites
  result.output = opts.output

proc validate_counts(cnt: counts, meta: SketchMeta, expected_autosomal_sites: int) =
  if cnt.sample_name != meta.sample_name:
    quit &"[somalier] sketch file '{meta.path}' changed sample name from '{meta.sample_name}' to '{cnt.sample_name}' while reading"
  if cnt.sites.len != expected_autosomal_sites:
    quit &"[somalier] sketch file '{meta.path}' has {cnt.sites.len} autosomal counts, expected {expected_autosomal_sites}"

proc call_charr_genotype(ac: allele_count, min_depth: int,
    hom_cutoff: float): int8 {.inline.} =
  let dp = int(ac.nref + ac.nalt)
  if dp < min_depth:
    return -1
  if ac.nother.float / (dp.float + ac.nother.float) > 0.04:
    return -1
  if ac.nref == 0'u32:
    return 2
  if ac.nalt == 0'u32:
    return 0
  let ab = ac.nalt.float / dp.float
  if ab <= hom_cutoff:
    return 0
  if ab >= (1.0 - hom_cutoff):
    return 2
  return -1

proc estimate_charr*(
    allele_counts: openArray[allele_count],
    pop_afs: openArray[float32],
    min_depth: int,
    hom_cutoff: float
  ): CharrStats =
  const min_contam_allele_af = 1e-6'f64
  doAssert allele_counts.len == pop_afs.len

  var total = 0'f64
  for i, ac in allele_counts:
    let gt = call_charr_genotype(ac, min_depth, hom_cutoff)
    let dp = ac.depth.float64
    let pB = pop_afs[i].float64
    var
      infiltrating = 0'f64
      contam_af = 0'f64

    case gt
    of 0'i8:
      result.n_hom_ref_usable.inc
      infiltrating = ac.nalt.float64
      contam_af = pB
    of 2'i8:
      result.n_hom_alt_usable.inc
      infiltrating = ac.nref.float64
      contam_af = 1'f64 - pB
    else:
      continue

    if contam_af <= min_contam_allele_af:
      if gt == 0'i8:
        result.n_hom_ref_usable.dec
      else:
        result.n_hom_alt_usable.dec
      continue

    total += infiltrating / (contam_af * dp)
    result.n_sites_usable.inc

  if result.n_sites_usable == 0:
    result.contamination = NaN
  else:
    result.contamination = total / result.n_sites_usable.float64

proc charr_rows(paths: seq[string], sites_path: string, min_depth: int,
    hom_cutoff: float = 0.15): seq[CharrRow] =
  let site_data = readSitesWithAF(sites_path)
  let sketches = read_sketch_headers(paths, site_data.pop_afs.len)
  var cnt: counts
  stderr.write_line &"[somalier] indexed {sketches.len} extracted files with {site_data.pop_afs.len} autosomal AF sites"

  for sketch in sketches:
    read_extracted(sketch.path, cnt)
    validate_counts(cnt, sketch, site_data.pop_afs.len)
    result.add((cnt.sample_name, estimate_charr(cnt.sites, site_data.pop_afs,
        min_depth, hom_cutoff)))

proc format_contamination(v: float64): string =
  if v.classify == fcNan:
    return "NaN"
  formatFloatClean(v.float32)

proc write_charr_rows(rows: openArray[CharrRow], output: string) =
  var fh: File
  let use_stdout = output in ["-", "stdout"]
  if use_stdout:
    fh = stdout
  elif not open(fh, output, fmWrite):
    quit "[somalier] couldn't open output file: " & output

  fh.write_line("#sample_name\tn_sites_usable\tcontamination_charr")
  for row in rows:
    fh.write_line(&"{row.sample_name}\t{row.stats.n_sites_usable}\t{format_contamination(row.stats.contamination)}")

  if not use_stdout:
    fh.close()

proc contamination_main*() =
  let cfg = parse_contamination_args(commandLineParams())
  let rows = charr_rows(cfg.extracted, cfg.sites, cfg.min_depth,
      cfg.charr_hom_cutoff)
  write_charr_rows(rows, cfg.output)
