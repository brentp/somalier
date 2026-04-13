import os
import streams
import strutils
import tables
import strformat
import math
import argparse
import ./common
import ./charr
import ./depthview
import ./pedrel

type ContaminationConfig = object
  sites: string
  min_depth: int
  charr_hom_rate: float
  charr_hom_alpha: float
  exact_hom_rate: float
  exact_hom_alpha: float
  output_prefix: string
  extracted: seq[string]

type PairContaminationStats* = object
  contamination*: float64
  n_sites_usable*: int

type CharrRow = tuple[sample_name: string, stats: CharrStats]
type PairRow = tuple[sample_name: string, anchor_sample: string, stats: PairContaminationStats]

type SketchMeta = object
  path: string
  sample_name: string
  n_sites: uint16
  nx_sites: uint16
  ny_sites: uint16

type ExactSiteData = object
  nref: int32
  nalt: int32
  usable: bool
  anchor_gt: int8

type LoadedSketch = object
  path: string
  sample_name: string
  sites: seq[allele_count]
  exact_sites: seq[ExactSiteData]

type PairSiteObservation = object
  receiver_nref: int
  receiver_nalt: int
  anchor_gt: int8
  pop_af: float64

type LoadedContaminationInputs = object
  pop_afs: seq[float32]
  sketches: seq[LoadedSketch]

const
  minLikelihoodProb = 1e-10'f64
  minPriorAf = 1e-6'f64
  defaultOutputPrefix = "somalier-contamination"
  exactGridStep = 0.01
  exactRefineTol = 1e-4
  exactMaxAlpha = 0.5
  defaultExactHomRate = 0.05
  defaultExactHomAlpha = 0.001
  goldenRatioConjugate = 0.3819660112501051

proc gt_prior*(pop_af: float64): array[3, float64]
proc site_log_likelihood*(nref, nalt: int, p: float64): float64
proc to_exact_site(ac: allele_count, min_depth: int,
    hom_minor_rate: float, hom_alpha: float): ExactSiteData

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
    # Contamination requires the exact same autosomal sketch ordering as the
    # sites file so both CHARR and the pairwise likelihood see aligned counts.
    if meta.n_sites.int != expected_autosomal_sites:
      quit &"[somalier] sketch file '{path}' has {meta.n_sites} autosomal sites, expected {expected_autosomal_sites} from --sites"
    result.add(meta)

proc normalize_output_prefix(output: string): string =
  result = if output.len == 0: defaultOutputPrefix else: output
  if result in ["-", "stdout"]:
    quit "[somalier] contamination writes multiple files; use -o PREFIX instead of stdout"
  if result.endsWith(".tsv"):
    result = result[0 ..< result.len - 4]
  if result.len == 0:
    result = defaultOutputPrefix

proc sample_output_path(prefix: string): string =
  prefix & ".samples.tsv"

proc pair_output_path(prefix: string): string =
  prefix & ".pairs.tsv"

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
    option("--charr-hom-rate", default = $defaultCharrHomRate,
      help = "tolerated minor-allele rate under the homozygous-like CHARR null model")
    option("--charr-hom-alpha", default = $defaultCharrHomAlpha,
      help = "minimum binomial tail probability required to keep a site as homozygous-like for CHARR")
    option("--exact-hom-rate", default = $defaultExactHomRate,
      help = "tolerated minor-allele rate under the stricter homozygous-like anchor model for exact contamination")
    option("--exact-hom-alpha", default = $defaultExactHomAlpha,
      help = "minimum binomial tail probability required to keep an anchor site as homozygous-like for exact contamination")
    option("-o", "--output", default = defaultOutputPrefix,
      help = "output prefix for PREFIX.samples.tsv and PREFIX.pairs.tsv")
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
    result.charr_hom_rate = parseFloat(opts.charr_hom_rate)
  except ValueError:
    quit "[somalier] --charr-hom-rate must be a number"
  if result.charr_hom_rate <= 0 or result.charr_hom_rate >= 0.5:
    quit "[somalier] --charr-hom-rate must be greater than 0 and less than 0.5"
  try:
    result.charr_hom_alpha = parseFloat(opts.charr_hom_alpha)
  except ValueError:
    quit "[somalier] --charr-hom-alpha must be a number"
  if result.charr_hom_alpha <= 0 or result.charr_hom_alpha >= 1:
    quit "[somalier] --charr-hom-alpha must be greater than 0 and less than 1"
  try:
    result.exact_hom_rate = parseFloat(opts.exact_hom_rate)
  except ValueError:
    quit "[somalier] --exact-hom-rate must be a number"
  if result.exact_hom_rate <= 0 or result.exact_hom_rate >= 0.5:
    quit "[somalier] --exact-hom-rate must be greater than 0 and less than 0.5"
  try:
    result.exact_hom_alpha = parseFloat(opts.exact_hom_alpha)
  except ValueError:
    quit "[somalier] --exact-hom-alpha must be a number"
  if result.exact_hom_alpha <= 0 or result.exact_hom_alpha >= 1:
    quit "[somalier] --exact-hom-alpha must be greater than 0 and less than 1"

  result.sites = opts.sites
  result.output_prefix = normalize_output_prefix(opts.output)

proc validate_counts(cnt: counts, meta: SketchMeta, expected_autosomal_sites: int) =
  if cnt.sample_name != meta.sample_name:
    quit &"[somalier] sketch file '{meta.path}' changed sample name from '{meta.sample_name}' to '{cnt.sample_name}' while reading"
  if cnt.sites.len != expected_autosomal_sites:
    quit &"[somalier] sketch file '{meta.path}' has {cnt.sites.len} autosomal counts, expected {expected_autosomal_sites}"

proc sample_site_usable(ac: allele_count, min_depth: int): bool {.inline.} =
  let dp = int(ac.nref + ac.nalt)
  if dp < min_depth:
    return false
  ac.nother.float / (dp.float + ac.nother.float) <= maxOtherFraction

proc load_contamination_inputs(paths: seq[string], sites_path: string,
    min_depth: int = 7, exact_hom_rate: float = defaultExactHomRate,
    exact_hom_alpha: float = defaultExactHomAlpha): LoadedContaminationInputs =
  let site_data = readSitesWithAF(sites_path)
  let headers = read_sketch_headers(paths, site_data.pop_afs.len)
  var cnt: counts
  result.pop_afs = site_data.pop_afs
  result.sketches = newSeqOfCap[LoadedSketch](headers.len)
  # The exact pairwise stage is all-vs-all, so loading each autosomal sketch
  # once avoids rereading the same files for every anchor pairing.
  for meta in headers:
    read_extracted(meta.path, cnt)
    validate_counts(cnt, meta, site_data.pop_afs.len)
    var exact_sites = newSeq[ExactSiteData](cnt.sites.len)
    for i, ac in cnt.sites:
      exact_sites[i] = to_exact_site(ac, min_depth, exact_hom_rate, exact_hom_alpha)
    result.sketches.add(LoadedSketch(
      path: meta.path,
      sample_name: cnt.sample_name,
      sites: cnt.sites,
      exact_sites: exact_sites
    ))
  stderr.write_line &"[somalier] loaded {result.sketches.len} extracted files with {result.pop_afs.len} autosomal AF sites"

proc gt_prior*(pop_af: float64): array[3, float64] =
  let q = clamp(pop_af, minPriorAf, 1'f64 - minPriorAf)
  [
    (1'f64 - q) * (1'f64 - q),
    2'f64 * q * (1'f64 - q),
    q * q
  ]

proc site_log_likelihood*(nref, nalt: int, p: float64): float64 =
  let pc = clamp(p, minLikelihoodProb, 1'f64 - minLikelihoodProb)
  nalt.float64 * ln(pc) + nref.float64 * ln(1'f64 - pc)

proc to_exact_site(ac: allele_count, min_depth: int,
    hom_minor_rate: float, hom_alpha: float): ExactSiteData =
  result.nref = ac.nref.int32
  result.nalt = ac.nalt.int32
  result.usable = sample_site_usable(ac, min_depth)
  if result.usable:
    var max_minor_by_depth = initTable[int, int]()
    result.anchor_gt = call_homozygous_like_genotype(ac, min_depth,
        hom_minor_rate, hom_alpha, max_minor_by_depth)
  else:
    result.anchor_gt = -1

proc logsumexp*(xs: openArray[float64]): float64 =
  if xs.len == 0:
    return NegInf
  var m = xs[0]
  for i in 1..<xs.len:
    if xs[i] > m:
      m = xs[i]
  if m.classify == fcNegInf:
    return m
  var s = 0'f64
  for x in xs:
    s += exp(x - m)
  m + ln(s)

proc anchor_expected_alt_prob*(anchor_gt: int8, alpha: float64, pop_af: float64): float64 =
  case anchor_gt
  of 0'i8:
    alpha * pop_af
  of 2'i8:
    1'f64 - alpha * (1'f64 - pop_af)
  else:
    NaN

proc pair_site_marginal_log_likelihood*(
    receiver_nref, receiver_nalt: int,
    anchor_gt: int8,
    alpha: float64,
    pop_af: float64
  ): float64 =
  let p = anchor_expected_alt_prob(anchor_gt, alpha, pop_af)
  if p.classify == fcNan:
    return NegInf
  site_log_likelihood(receiver_nref, receiver_nalt, p)

proc pair_total_log_likelihood*(sites: openArray[PairSiteObservation], alpha: float64): float64 =
  for site in sites:
    result += pair_site_marginal_log_likelihood(site.receiver_nref, site.receiver_nalt,
        site.anchor_gt, alpha, site.pop_af)

proc pair_total_log_likelihood(receiver, anchor: LoadedSketch,
    pop_afs: openArray[float32], alpha: float64): float64 =
  for i in 0..<receiver.exact_sites.len:
    let rs = receiver.exact_sites[i]
    let `as` = anchor.exact_sites[i]
    if not (rs.usable and `as`.usable):
      continue
    if `as`.anchor_gt notin [0'i8, 2'i8]:
      continue
    result += pair_site_marginal_log_likelihood(rs.nref.int, rs.nalt.int,
        `as`.anchor_gt, alpha, pop_afs[i].float64)

proc golden_section_max(receiver, anchor: LoadedSketch,
    pop_afs: openArray[float32],
    lo_in, hi_in: float64): tuple[alpha: float64, ll: float64] =
  var
    lo = lo_in
    hi = hi_in
    c = hi - goldenRatioConjugate * (hi - lo)
    d = lo + goldenRatioConjugate * (hi - lo)
    fc = pair_total_log_likelihood(receiver, anchor, pop_afs, c)
    fd = pair_total_log_likelihood(receiver, anchor, pop_afs, d)

  while hi - lo > exactRefineTol:
    if fc < fd:
      lo = c
      c = d
      fc = fd
      d = lo + goldenRatioConjugate * (hi - lo)
      fd = pair_total_log_likelihood(receiver, anchor, pop_afs, d)
    else:
      hi = d
      d = c
      fd = fc
      c = hi - goldenRatioConjugate * (hi - lo)
      fc = pair_total_log_likelihood(receiver, anchor, pop_afs, c)

  if fc >= fd:
    (c, fc)
  else:
    (d, fd)

proc golden_section_max(sites: openArray[PairSiteObservation],
    lo_in, hi_in: float64): tuple[alpha: float64, ll: float64] =
  var
    lo = lo_in
    hi = hi_in
    c = hi - goldenRatioConjugate * (hi - lo)
    d = lo + goldenRatioConjugate * (hi - lo)
    fc = pair_total_log_likelihood(sites, c)
    fd = pair_total_log_likelihood(sites, d)

  while hi - lo > exactRefineTol:
    if fc < fd:
      lo = c
      c = d
      fc = fd
      d = lo + goldenRatioConjugate * (hi - lo)
      fd = pair_total_log_likelihood(sites, d)
    else:
      hi = d
      d = c
      fd = fc
      c = hi - goldenRatioConjugate * (hi - lo)
      fc = pair_total_log_likelihood(sites, c)

  if fc >= fd:
    (c, fc)
  else:
    (d, fd)

proc estimate_pair_contamination*(
    receiver_counts: openArray[allele_count],
    anchor_counts: openArray[allele_count],
    pop_afs: openArray[float32],
    min_depth: int,
    exact_hom_rate: float = defaultExactHomRate,
    exact_hom_alpha: float = defaultExactHomAlpha
  ): PairContaminationStats =
  doAssert receiver_counts.len == anchor_counts.len
  doAssert receiver_counts.len == pop_afs.len

  var usable_sites = newSeqOfCap[PairSiteObservation](receiver_counts.len)
  var max_minor_by_depth = initTable[int, int]()
  for i in 0..<receiver_counts.len:
    if not (sample_site_usable(receiver_counts[i], min_depth) and
        sample_site_usable(anchor_counts[i], min_depth)):
      continue
    let anchor_gt = call_homozygous_like_genotype(anchor_counts[i], min_depth,
        exact_hom_rate, exact_hom_alpha, max_minor_by_depth)
    if anchor_gt notin [0'i8, 2'i8]:
      continue
    usable_sites.add(PairSiteObservation(
      receiver_nref: receiver_counts[i].nref.int,
      receiver_nalt: receiver_counts[i].nalt.int,
      anchor_gt: anchor_gt,
      pop_af: pop_afs[i].float64
    ))

  result.n_sites_usable = usable_sites.len
  if usable_sites.len == 0:
    result.contamination = NaN
    return

  var
    best_alpha = 0'f64
    best_ll = NegInf
  let n_steps = int(exactMaxAlpha / exactGridStep)
  for step in 0..n_steps:
    let alpha = min(exactMaxAlpha, step.float64 * exactGridStep)
    let ll = pair_total_log_likelihood(usable_sites, alpha)
    if ll > best_ll:
      best_ll = ll
      best_alpha = alpha

  let refine_lo = max(0'f64, best_alpha - exactGridStep)
  let refine_hi = min(exactMaxAlpha, best_alpha + exactGridStep)
  if refine_hi - refine_lo <= exactRefineTol:
    result.contamination = best_alpha
    return

  let refined = golden_section_max(usable_sites, refine_lo, refine_hi)
  if refined.ll > best_ll:
    result.contamination = refined.alpha
  else:
    result.contamination = best_alpha

proc estimate_pair_contamination(receiver, anchor: LoadedSketch,
    pop_afs: openArray[float32]): PairContaminationStats =
  doAssert receiver.exact_sites.len == anchor.exact_sites.len
  doAssert receiver.exact_sites.len == pop_afs.len

  for i in 0..<receiver.exact_sites.len:
    if receiver.exact_sites[i].usable and anchor.exact_sites[i].usable and
        anchor.exact_sites[i].anchor_gt in [0'i8, 2'i8]:
      result.n_sites_usable.inc

  if result.n_sites_usable == 0:
    result.contamination = NaN
    return

  var
    best_alpha = 0'f64
    best_ll = NegInf
  let n_steps = int(exactMaxAlpha / exactGridStep)
  for step in 0..n_steps:
    let alpha = min(exactMaxAlpha, step.float64 * exactGridStep)
    let ll = pair_total_log_likelihood(receiver, anchor, pop_afs, alpha)
    if ll > best_ll:
      best_ll = ll
      best_alpha = alpha

  let refine_lo = max(0'f64, best_alpha - exactGridStep)
  let refine_hi = min(exactMaxAlpha, best_alpha + exactGridStep)
  if refine_hi - refine_lo <= exactRefineTol:
    result.contamination = best_alpha
    return

  let refined = golden_section_max(receiver, anchor, pop_afs, refine_lo, refine_hi)
  if refined.ll > best_ll:
    result.contamination = refined.alpha
  else:
    result.contamination = best_alpha

proc charr_rows(sketches: openArray[LoadedSketch], pop_afs: openArray[float32],
    min_depth: int, hom_minor_rate: float, hom_alpha: float): seq[CharrRow] =
  for sketch in sketches:
    result.add((sketch.sample_name, estimate_charr(sketch.sites, pop_afs,
        min_depth, hom_minor_rate, hom_alpha)))

proc pair_rows(sketches: openArray[LoadedSketch],
    pop_afs: openArray[float32]): seq[PairRow] =
  for i, receiver in sketches:
    # Pair rows are directional, so we hold the receiver fixed and evaluate
    # every other sample as the candidate matched anchor.
    for j, anchor in sketches:
      if i == j:
        continue
      result.add((receiver.sample_name, anchor.sample_name,
          estimate_pair_contamination(receiver, anchor, pop_afs)))

proc format_contamination(v: float64): string =
  if v.classify == fcNan:
    return "NaN"
  formatFloatClean(v.float32)

proc write_charr_rows(rows: openArray[CharrRow], output: string) =
  var fh: File
  if not open(fh, output, fmWrite):
    quit "[somalier] couldn't open output file: " & output
  defer:
    fh.close()

  fh.write_line("#sample_name\tn_sites_usable\tcontamination_charr")
  for row in rows:
    fh.write_line(&"{row.sample_name}\t{row.stats.n_sites_usable}\t{format_contamination(row.stats.contamination)}")

proc write_pair_rows(rows: openArray[PairRow], output: string) =
  var fh: File
  if not open(fh, output, fmWrite):
    quit "[somalier] couldn't open output file: " & output
  defer:
    fh.close()

  fh.write_line("#sample_name\tanchor_sample\tn_sites_usable\tcontamination_mle")
  for row in rows:
    fh.write_line(&"{row.sample_name}\t{row.anchor_sample}\t{row.stats.n_sites_usable}\t{format_contamination(row.stats.contamination)}")

proc contamination_main*() =
  let cfg = parse_contamination_args(commandLineParams())
  let loaded = load_contamination_inputs(cfg.extracted, cfg.sites, cfg.min_depth,
      cfg.exact_hom_rate, cfg.exact_hom_alpha)
  let sample_rows = charr_rows(loaded.sketches, loaded.pop_afs, cfg.min_depth,
      cfg.charr_hom_rate, cfg.charr_hom_alpha)
  let pairs = pair_rows(loaded.sketches, loaded.pop_afs)
  let sample_path = sample_output_path(cfg.output_prefix)
  let pair_path = pair_output_path(cfg.output_prefix)
  write_charr_rows(sample_rows, sample_path)
  write_pair_rows(pairs, pair_path)
  stderr.write_line &"[somalier] wrote {sample_path} and {pair_path}"
