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
  exact_error_rate: float
  output_prefix: string
  extracted: seq[string]
  tumor_normal_pair_mode: bool

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
  log_contam_priors: array[3, float64]

type GridSearchResult = object
  alpha: float64
  ll: float64
  lo: float64
  hi: float64

type LoadedContaminationInputs = object
  pop_afs: seq[float32]
  log_priors: seq[array[3, float64]]
  sketches: seq[LoadedSketch]

const
  # Clamp site probabilities away from exact 0/1 before taking logs.
  minLikelihoodProb = 1e-10'f64
  # Keep Hardy-Weinberg genotype priors finite at rare/common AF extremes.
  minPriorAf = 1e-6'f64
  # Shared fallback prefix when the caller omits -o/--output.
  defaultOutputPrefix = "somalier-contamination"
  # Stop refining once the local search bracket is this narrow.
  exactRefineTol = 1e-4
  # Search contamination fractions on [0, exactMaxAlpha] for Conpair-like cases
  # that can exceed 50% apparent contamination in heavily distorted tumors.
  exactMaxAlpha = 1.0
  # Stricter homozygous-like anchor filter than CHARR uses.
  defaultExactHomRate = 0.05
  # Tail-probability cutoff for keeping an anchor site homozygous-like.
  defaultExactHomAlpha = 0.001
  # Approximate per-read sequencing error used in the exact likelihood.
  defaultExactErrorRate = 0.002
  # Constant used by golden-section search to shrink the refinement interval.
  goldenRatioConjugate = 0.3819660112501051

proc site_log_likelihood*(nref, nalt: int, p: float64): float64
proc to_exact_site(ac: allele_count, min_depth: int,
    hom_minor_rate: float, hom_alpha: float,
    max_minor_by_depth: var Table[int, int]): ExactSiteData

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
  # Read just the fixed header so we can validate many sketches cheaply first.
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
    # Accept either a bare prefix or a legacy-looking .tsv path from callers.
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
    option("--min-depth", default = "15",
      help = "minimum depth required for a site to be usable")
    option("--charr-hom-rate", default = $defaultCharrHomRate,
      help = "tolerated minor-allele rate under the homozygous-like CHARR null model")
    option("--charr-hom-alpha", default = $defaultCharrHomAlpha,
      help = "minimum binomial tail probability required to keep a site as homozygous-like for CHARR")
    option("--exact-hom-rate", default = $defaultExactHomRate,
      help = "tolerated minor-allele rate under the stricter homozygous-like anchor model for exact contamination")
    option("--exact-hom-alpha", default = $defaultExactHomAlpha,
      help = "minimum binomial tail probability required to keep an anchor site as homozygous-like for exact contamination")
    option("--exact-error-rate", default = $defaultExactErrorRate,
      help = "per-read sequencing error rate used by the exact pairwise likelihood")
    flag("-p", "--tumor-normal-pair",
      help = "emit a single tumor-target normal-reference contamination line using exactly 2 extracted .somalier inputs; the first is assumed to be the tumor")
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
  result.tumor_normal_pair_mode = opts.tumor_normal_pair
  if result.tumor_normal_pair_mode:
    if result.extracted.len != 2:
      quit "[somalier] --tumor-normal-pair requires exactly 2 extracted .somalier files; the first is the tumor and the second is the normal"
    for path in result.extracted:
      if not fileExists(path):
        quit "[somalier] unable to open extracted somalier file: " & path
      if not path.endsWith(".somalier"):
        quit "[somalier] --tumor-normal-pair requires exactly 2 extracted .somalier files"
  else:
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
  try:
    result.exact_error_rate = parseFloat(opts.exact_error_rate)
  except ValueError:
    quit "[somalier] --exact-error-rate must be a number"
  if result.exact_error_rate < 0 or result.exact_error_rate >= 0.5:
    quit "[somalier] --exact-error-rate must be at least 0 and less than 0.5"

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

proc gt_prior*(pop_af: float64): array[3, float64] =
  let q = clamp(pop_af, minPriorAf, 1'f64 - minPriorAf)
  # Use Hardy-Weinberg priors so contaminant genotypes follow the site AF.
  [
    (1'f64 - q) * (1'f64 - q),
    2'f64 * q * (1'f64 - q),
    q * q
  ]

proc log_gt_prior(pop_af: float64): array[3, float64] =
  let priors = gt_prior(pop_af)
  [
    ln(priors[0]),
    ln(priors[1]),
    ln(priors[2])
  ]

proc compute_log_gt_priors(pop_afs: openArray[float32]): seq[array[3, float64]] =
  result = newSeq[array[3, float64]](pop_afs.len)
  for i, pop_af in pop_afs:
    result[i] = log_gt_prior(pop_af.float64)

proc load_contamination_inputs(paths: seq[string], sites_path: string,
    min_depth: int = 15, exact_hom_rate: float = defaultExactHomRate,
    exact_hom_alpha: float = defaultExactHomAlpha): LoadedContaminationInputs =
  let site_data = readSitesWithAF(sites_path)
  let headers = read_sketch_headers(paths, site_data.pop_afs.len)
  var cnt: counts
  result.pop_afs = site_data.pop_afs
  result.log_priors = compute_log_gt_priors(result.pop_afs)
  result.sketches = newSeqOfCap[LoadedSketch](headers.len)
  # The exact pairwise stage is all-vs-all, so loading each autosomal sketch
  # once avoids rereading the same files for every anchor pairing.
  for meta in headers:
    read_extracted(meta.path, cnt)
    validate_counts(cnt, meta, site_data.pop_afs.len)
    var exact_sites = newSeq[ExactSiteData](cnt.sites.len)
    var max_minor_by_depth = initTable[int, int]()
    for i, ac in cnt.sites:
      # Precompute exact-stage site filters once so every anchor pairing can reuse them.
      exact_sites[i] = to_exact_site(ac, min_depth, exact_hom_rate,
          exact_hom_alpha, max_minor_by_depth)
    max_minor_by_depth.clear()
    result.sketches.add(LoadedSketch(
      path: meta.path,
      sample_name: cnt.sample_name,
      sites: cnt.sites,
      exact_sites: exact_sites
    ))
  stderr.write_line &"[somalier] loaded {result.sketches.len} extracted files with {result.pop_afs.len} autosomal AF sites"

proc site_log_likelihood*(nref, nalt: int, p: float64): float64 =
  let pc = clamp(p, minLikelihoodProb, 1'f64 - minLikelihoodProb)
  nalt.float64 * ln(pc) + nref.float64 * ln(1'f64 - pc)

proc to_exact_site(ac: allele_count, min_depth: int,
    hom_minor_rate: float, hom_alpha: float,
    max_minor_by_depth: var Table[int, int]): ExactSiteData =
  result.nref = ac.nref.int32
  result.nalt = ac.nalt.int32
  result.usable = sample_site_usable(ac, min_depth)
  result.anchor_gt = if result.usable:
    # Cache the anchor's homozygous-like status so pairwise scoring stays lightweight.
    call_homozygous_like_genotype(ac, min_depth,
        hom_minor_rate, hom_alpha, max_minor_by_depth)
  else:
    -1'i8

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
  # Sum in log-space so mixed genotype likelihoods stay numerically stable.
  m + ln(s)

proc genotype_alt_fraction(gt: int8): float64 {.inline.} =
  case gt
  of 0'i8:
    0'f64
  of 1'i8:
    0.5'f64
  of 2'i8:
    1'f64
  else:
    NaN

proc observed_alt_prob(clean_gt, contam_gt: int8, alpha: float64,
    error_rate: float64): float64 =
  let clean_alt = genotype_alt_fraction(clean_gt)
  let contam_alt = genotype_alt_fraction(contam_gt)
  if clean_alt.classify == fcNan or contam_alt.classify == fcNan:
    return NaN
  # Mix clean and contaminant allele fractions before folding in sequencing error.
  let latent_alt = (1'f64 - alpha) * clean_alt + alpha * contam_alt
  latent_alt * (1'f64 - error_rate) + (1'f64 - latent_alt) * error_rate

proc pair_site_marginal_log_likelihood*(
    receiver_nref, receiver_nalt: int,
    anchor_gt: int8,
    alpha: float64,
    log_contam_priors: array[3, float64],
    error_rate: float64
  ): float64 =
  if anchor_gt notin [0'i8, 2'i8]:
    return NegInf
  var terms: array[3, float64]
  for contam_gt in 0'i8..2'i8:
    # Marginalize over contaminant genotype instead of collapsing to mean AF first.
    let p = observed_alt_prob(anchor_gt, contam_gt, alpha, error_rate)
    terms[contam_gt.int] = log_contam_priors[contam_gt.int] +
        site_log_likelihood(receiver_nref, receiver_nalt, p)
  logsumexp(terms)

proc pair_total_log_likelihood*(sites: openArray[PairSiteObservation], alpha: float64,
    error_rate: float64): float64 =
  for site in sites:
    result += pair_site_marginal_log_likelihood(site.receiver_nref, site.receiver_nalt,
        site.anchor_gt, alpha, site.log_contam_priors, error_rate)

proc build_grid(lo_in, hi_in, step: float64): seq[float64] =
  let
    lo = max(0'f64, min(exactMaxAlpha, lo_in))
    hi = max(0'f64, min(exactMaxAlpha, hi_in))
  if hi < lo or step <= 0'f64:
    return @[]
  var alpha = lo
  while alpha < hi:
    result.add(alpha)
    alpha += step
  if result.len == 0 or result[^1] != hi:
    result.add(hi)

proc search_grid(sites: openArray[PairSiteObservation], error_rate: float64,
    values: openArray[float64]): GridSearchResult =
  doAssert values.len > 0
  for i in 1..<values.len:
    doAssert values[i] >= values[i - 1]

  var best_idx = 0
  result.alpha = values[0]
  result.ll = pair_total_log_likelihood(sites, values[0], error_rate)
  for i in 1..<values.len:
    let ll = pair_total_log_likelihood(sites, values[i], error_rate)
    if ll > result.ll:
      result.ll = ll
      result.alpha = values[i]
      best_idx = i

  result.lo = if best_idx == 0: values[0] else: values[best_idx - 1]
  result.hi = if best_idx == values.high: values[^1] else: values[best_idx + 1]

proc maybe_take_better(best: var GridSearchResult, candidate: GridSearchResult) =
  if candidate.ll > best.ll:
    best = candidate
    return
  if candidate.ll == best.ll and candidate.alpha == best.alpha and
      (candidate.hi - candidate.lo) < (best.hi - best.lo):
    best = candidate

proc collect_pair_sites(receiver, anchor: LoadedSketch,
    log_priors: openArray[array[3, float64]]): seq[PairSiteObservation] =
  for i in 0..<receiver.exact_sites.len:
    let rs = receiver.exact_sites[i]
    let `as` = anchor.exact_sites[i]
    if not (rs.usable and `as`.usable):
      continue
    if `as`.anchor_gt notin [0'i8, 2'i8]:
      continue
    # Reuse per-site priors computed at load time so pair scoring only varies alpha.
    result.add(PairSiteObservation(
      receiver_nref: rs.nref.int,
      receiver_nalt: rs.nalt.int,
      anchor_gt: `as`.anchor_gt,
      log_contam_priors: log_priors[i]
    ))

proc golden_section_max(sites: openArray[PairSiteObservation],
    error_rate: float64,
    lo_in, hi_in: float64): tuple[alpha: float64, ll: float64] =
  var
    lo = lo_in
    hi = hi_in
    c = hi - goldenRatioConjugate * (hi - lo)
    d = lo + goldenRatioConjugate * (hi - lo)
    fc = pair_total_log_likelihood(sites, c, error_rate)
    fd = pair_total_log_likelihood(sites, d, error_rate)

  # Use the same refinement path for the standalone test helper variant.
  while hi - lo > exactRefineTol:
    if fc < fd:
      lo = c
      c = d
      fc = fd
      d = lo + goldenRatioConjugate * (hi - lo)
      fd = pair_total_log_likelihood(sites, d, error_rate)
    else:
      hi = d
      d = c
      fd = fc
      c = hi - goldenRatioConjugate * (hi - lo)
      fc = pair_total_log_likelihood(sites, c, error_rate)

  if fc >= fd:
    (c, fc)
  else:
    (d, fd)

proc estimate_pair_contamination_from_sites(
    usable_sites: openArray[PairSiteObservation],
    exact_error_rate: float64
  ): PairContaminationStats =
  result.n_sites_usable = usable_sites.len
  if usable_sites.len == 0:
    result.contamination = NaN
    return

  var best = search_grid(usable_sites, exact_error_rate, [
    0'f64, 0.01'f64, 0.02'f64, 0.05'f64, 0.1'f64,
    0.25'f64, 0.5'f64, 0.75'f64, 1.0'f64
  ])

  if best.alpha >= 0.75'f64:
    let high_coarse = search_grid(usable_sites, exact_error_rate,
        build_grid(0.75'f64, exactMaxAlpha, 0.05'f64))
    maybe_take_better(best, high_coarse)
    let high_fine = search_grid(usable_sites, exact_error_rate,
        build_grid(best.lo, best.hi, 0.01'f64))
    maybe_take_better(best, high_fine)

  if best.hi - best.lo <= exactRefineTol:
    result.contamination = best.alpha
    return

  let refined = golden_section_max(usable_sites, exact_error_rate, best.lo, best.hi)
  if refined.ll > best.ll:
    result.contamination = refined.alpha
  else:
    result.contamination = best.alpha

proc estimate_pair_contamination*(
    receiver_counts: openArray[allele_count],
    anchor_counts: openArray[allele_count],
    pop_afs: openArray[float32],
    min_depth: int,
    exact_hom_rate: float = defaultExactHomRate,
    exact_hom_alpha: float = defaultExactHomAlpha,
    exact_error_rate: float = defaultExactErrorRate
  ): PairContaminationStats =
  doAssert receiver_counts.len == anchor_counts.len
  doAssert receiver_counts.len == pop_afs.len

  var usable_sites = newSeqOfCap[PairSiteObservation](receiver_counts.len)
  var max_minor_by_depth = initTable[int, int]()
  let log_priors = compute_log_gt_priors(pop_afs)
  for i in 0..<receiver_counts.len:
    if not (sample_site_usable(receiver_counts[i], min_depth) and
        sample_site_usable(anchor_counts[i], min_depth)):
      continue
    # Anchor-side gating removes sites that are already too imbalanced to anchor on.
    let anchor_gt = call_homozygous_like_genotype(anchor_counts[i], min_depth,
        exact_hom_rate, exact_hom_alpha, max_minor_by_depth)
    if anchor_gt notin [0'i8, 2'i8]:
      continue
    usable_sites.add(PairSiteObservation(
      receiver_nref: receiver_counts[i].nref.int,
      receiver_nalt: receiver_counts[i].nalt.int,
      anchor_gt: anchor_gt,
      log_contam_priors: log_priors[i]
    ))
  result = estimate_pair_contamination_from_sites(usable_sites, exact_error_rate)

proc estimate_pair_contamination(receiver, anchor: LoadedSketch,
    log_priors: openArray[array[3, float64]],
    exact_error_rate: float): PairContaminationStats =
  doAssert receiver.exact_sites.len == anchor.exact_sites.len
  doAssert receiver.exact_sites.len == log_priors.len

  let usable_sites = collect_pair_sites(receiver, anchor, log_priors)
  # Reuse the precomputed per-sample site metadata for the full all-vs-all path.
  result = estimate_pair_contamination_from_sites(usable_sites, exact_error_rate)

proc charr_rows(sketches: openArray[LoadedSketch], pop_afs: openArray[float32],
    min_depth: int, hom_minor_rate: float, hom_alpha: float): seq[CharrRow] =
  for sketch in sketches:
    # CHARR remains sample-level, so each sketch is scored independently here.
    result.add((sketch.sample_name, estimate_charr(sketch.sites, pop_afs,
        min_depth, hom_minor_rate, hom_alpha)))

proc pair_rows(sketches: openArray[LoadedSketch],
    log_priors: openArray[array[3, float64]],
    exact_error_rate: float): seq[PairRow] =
  for i, receiver in sketches:
    # Pair rows are directional, so we hold the receiver fixed and evaluate
    # every other sample as the candidate matched anchor.
    for j, anchor in sketches:
      if i == j:
        continue
      result.add((receiver.sample_name, anchor.sample_name,
          estimate_pair_contamination(receiver, anchor, log_priors, exact_error_rate)))

proc format_contamination(v: float64): string =
  if v.classify == fcNan:
    return "NaN"
  # Keep TSV output compact and stable with the same float formatting as relate.
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

proc pair_row_line(row: PairRow): string =
  &"{row.sample_name}\t{row.anchor_sample}\t{row.stats.n_sites_usable}\t{format_contamination(row.stats.contamination)}"

proc tumor_normal_pair_row(loaded: LoadedContaminationInputs,
    exact_error_rate: float): PairRow =
  doAssert loaded.sketches.len == 2
  let stats = estimate_pair_contamination(loaded.sketches[0], loaded.sketches[1],
      loaded.log_priors, exact_error_rate)
  (loaded.sketches[0].sample_name, loaded.sketches[1].sample_name, stats)

proc contamination_main*() =
  let cfg = parse_contamination_args(commandLineParams())
  let loaded = load_contamination_inputs(cfg.extracted, cfg.sites, cfg.min_depth,
      cfg.exact_hom_rate, cfg.exact_hom_alpha)
  if cfg.tumor_normal_pair_mode:
    if cfg.output_prefix != defaultOutputPrefix:
      quit "[somalier] --tumor-normal-pair writes a single line to stdout; omit -o/--output"
    let row = tumor_normal_pair_row(loaded, cfg.exact_error_rate)
    stdout.write_line(pair_row_line(row))
    return
  let sample_rows = charr_rows(loaded.sketches, loaded.pop_afs, cfg.min_depth,
      cfg.charr_hom_rate, cfg.charr_hom_alpha)
  let pairs = pair_rows(loaded.sketches, loaded.log_priors, cfg.exact_error_rate)
  let sample_path = sample_output_path(cfg.output_prefix)
  let pair_path = pair_output_path(cfg.output_prefix)
  write_charr_rows(sample_rows, sample_path)
  write_pair_rows(pairs, pair_path)
  stderr.write_line &"[somalier] wrote {sample_path} and {pair_path}"
