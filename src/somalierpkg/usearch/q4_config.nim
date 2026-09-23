## Compile-time configuration for the Q4 HNSW relatedness prefilter.
##
## Defaults are the validated operating point. Advanced builds may
## override values with `-d:name=value`. Print every value for each run and,
## when persistent indexes are added, store them in metadata so incompatible
## builds cannot mix indexes.
##
## Q4/HNSW execution is compiled into Somalier only with
## `-d:somalier_usearch`. This module remains dependency-free so an
## exhaustive-only build can still resolve the automatic-mode boundary and
## issue a useful rebuild error without importing or linking USearch.

const
  ## Pair-generation mode compiled into `somalier relate`. `auto` uses
  ## exhaustive scoring through `exhaustive_max_samples` and Q4/HNSW above
  ## that limit. `exhaustive` and `q4-hnsw` force one implementation for every
  ## cohort size; forcing Q4 is useful for validation on small cohorts.
  ## Override with `-d:somalier_q4_candidate_mode=q4-hnsw` or
  ## `-d:somalier_q4_candidate_mode=exhaustive`.
  candidate_mode* {.strdefine: "somalier_q4_candidate_mode".} = "auto"

  ## Largest cohort scored exhaustively when `candidate_mode` is `auto`.
  ## Q4/HNSW is selected only when the number of expanded `.somalier` inputs
  ## is greater than this value, so the default keeps exactly 5,000 samples on
  ## the exhaustive path and switches at 5,001. Override with
  ## `-d:somalier_q4_exhaustive_max_samples=10000`.
  exhaustive_max_samples* {.
    intdefine: "somalier_q4_exhaustive_max_samples".} = 5000

  ## Number of samples projected in one matrix-multiplication batch. This
  ## bounds transient standardized-genotype and projection memory without
  ## changing sketches or candidate semantics. Larger batches may improve
  ## throughput and consume proportionally more working memory. Override with
  ## `-d:somalier_q4_projection_batch_size=512`.
  projection_batch_size* {.
    intdefine: "somalier_q4_projection_batch_size".} = 256

  ## Number of concurrent USearch queries and the search-concurrency capacity
  ## reserved by the index. It is capped at the sample count at runtime. This
  ## affects wall time and temporary query state, not candidate semantics.
  ## Override with `-d:somalier_q4_threads=4`.
  query_threads* {.intdefine: "somalier_q4_threads".} = 1

  ## Largest Q4/HNSW cohort for which Somalier writes interactive HTML. A
  ## value of zero disables sparse-mode HTML. The TSV outputs are unaffected.
  ## Override with `-d:somalier_q4_html_max_samples=0`.
  html_max_samples* {.intdefine: "somalier_q4_html_max_samples".} = 10_000

  ## Optional path receiving every Q4-admitted pair, its packed-Q4 cosine, and
  ## `direct` or `rescue` admission reason. The empty default disables this
  ## diagnostic output. It may be very large and is intended for validation.
  ## Override with `-d:somalier_q4_candidate_output=/tmp/q4-candidates.tsv`.
  candidate_output_path* {.
    strdefine: "somalier_q4_candidate_output".} = ""

  ## Number of coordinates in each random projection. This is also the
  ## compile-time dimension used by `usearch_q4`, so its define intentionally
  ## has the wrapper's name rather than a Somalier-specific name. It must be
  ## positive and even because two signed Q4 coordinates are packed per byte.
  ## A record occupies `dimensions / 2 + 4` bytes: 4,096 dimensions use 2,052
  ## bytes including the stored squared norm. More dimensions increase index
  ## size and projection work linearly and may improve neighbor recall.
  ## Override with `-d:usearchQ4Dimensions=3072`.
  projection_dimensions* {.intdefine: "usearchQ4Dimensions".} = 4096

  ## Seed for the row-major `splitmix64-rademacher-v1` projection matrix. The
  ## same ordered panel, frequencies, dimensions, and seed deterministically
  ## produce the same projection. Changing the seed changes every Q4 sketch
  ## and therefore requires rebuilding the index. Values must be nonnegative.
  ## Override with `-d:somalier_q4_projection_seed=2718`.
  projection_seed_int* {.intdefine: "somalier_q4_projection_seed".} = 1729

  ## HNSW graph connectivity (`M`). Larger values retain more graph links per
  ## sample and generally improve search recall at the cost of index memory,
  ## construction time, and some query work. Override with
  ## `-d:somalier_q4_m=48`.
  hnsw_connectivity* {.intdefine: "somalier_q4_m".} = 32

  ## Size of the dynamic candidate list used while inserting an item into the
  ## HNSW graph (`efConstruction`). Larger values can build a higher-quality
  ## graph but increase construction time and temporary work. They do not
  ## change the packed Q4 vector size. Override with
  ## `-d:somalier_q4_ef_construction=500`.
  hnsw_ef_construction* {.intdefine: "somalier_q4_ef_construction".} = 400

  ## Size of the HNSW search expansion used for each top-K query (`efSearch`).
  ## Larger values usually improve approximate-neighbor recall and increase
  ## query time. This should normally be at least `neighbor_count`. Override
  ## with `-d:somalier_q4_ef_search=120`.
  hnsw_ef_search* {.intdefine: "somalier_q4_ef_search".} = 80

  ## Number of non-self neighbors requested for each sample. With reciprocal
  ## admission, the undirected graph contains at most `N * topK / 2` pairs;
  ## union admission contains at most `N * topK`. Raising top-K can recover
  ## relatives hidden by dense local neighborhoods, while increasing query
  ## output, neighbor-list memory, and exact candidate checks. Override with
  ## `-d:somalier_q4_top_k=60`.
  neighbor_count* {.intdefine: "somalier_q4_top_k".} = 40

  ## Candidate queue rule, encoded as an integer compile-time define. `1`
  ## admits an edge only when both samples retrieve each other; `0` admits the
  ## union of directed top-K results. Reciprocal mode is substantially more
  ## selective. Override with `-d:somalier_q4_reciprocal=0` for union mode.
  require_reciprocal_int* {.intdefine: "somalier_q4_reciprocal".} = 1

  ## Exact packed-Q4 cosine at or above this value is sent directly to exact
  ## Somalier pair scoring. The value is in thousandths because Nim's numeric
  ## compile-time defines are integers: `200` means cosine 0.200. Lowering it
  ## usually increases recall and candidate count. Override with
  ## `-d:somalier_q4_direct_cosine_milli=180`.
  direct_cosine_milli* {.intdefine: "somalier_q4_direct_cosine_milli".} = 200

  ## Lower edge of the positional-rescue band, in cosine thousandths. Pairs
  ## below this value are discarded. Pairs from this value up to, but not
  ## including, `direct_cosine_milli` must pass the positional rescue test. It
  ## must be lower than the direct threshold. Override with
  ## `-d:somalier_q4_rescue_cosine_milli=90`.
  rescue_cosine_milli* {.intdefine: "somalier_q4_rescue_cosine_milli".} = 110

  ## Absolute clipping limit applied after a projected vector is normalized
  ## and multiplied by `sqrt(dimensions)`. The value is in thousandths, so
  ## `3500` means 3.5. A lower limit saturates more coordinates; a higher limit
  ## preserves more tails but leaves fewer Q4 levels near zero when scale is
  ## reduced to fit. Override with `-d:somalier_q4_clip_milli=3000`.
  quantization_clip_milli* {.intdefine: "somalier_q4_clip_milli".} = 3500

  ## Multiplier applied after clipping and before rounding to an integer Q4
  ## coordinate. The value is in thousandths, so `2000` means multiply by 2.
  ## `clip * scale` may not exceed 7, the largest signed Q4 magnitude used by
  ## this format. Changing clip or scale changes cosine geometry and requires
  ## rebuilding and revalidating the index. Override with
  ## `-d:somalier_q4_scale_milli=1750`.
  quantization_scale_milli* {.intdefine: "somalier_q4_scale_milli".} = 2000

  ## Strict upper allele-balance boundary for a homozygous-reference call, in
  ## thousandths. The default `10` means `AB < 0.01`. Production Q4 projection,
  ## rescue, and exact scoring must all consume bitsets made with this boundary.
  ## Override with `-d:somalier_q4_hom_ref_cutoff_milli=20`.
  hom_ref_cutoff_milli* {.intdefine: "somalier_q4_hom_ref_cutoff_milli".} = 10

  ## Strict lower allele-balance boundary for a homozygous-alternate call, in
  ## thousandths. The default `990` means `AB > 0.99`. Values between the two
  ## homozygous boundaries are classified using `--min-ab` or remain unknown.
  ## Override with `-d:somalier_q4_hom_alt_cutoff_milli=980`.
  hom_alt_cutoff_milli* {.intdefine: "somalier_q4_hom_alt_cutoff_milli".} = 990

  ## Number of consecutive ordered panel sites in each positional-rescue
  ## window. Windows never cross chromosomes. Each of the two tilings advances
  ## by this width. Wider windows are more specific but are more vulnerable to
  ## errors and missing calls. Support counts below are derived from this width
  ## and their permille settings. Override with
  ## `-d:somalier_q4_rescue_width=128`.
  rescue_width* {.intdefine: "somalier_q4_rescue_width".} = 176

  ## Starting-site offset of the second rescue-window tiling relative to the
  ## first tiling at offset zero. It must be from zero through width minus one.
  ## Half the width covers boundaries between windows in the first tiling.
  ## Override with `-d:somalier_q4_rescue_offset=64` when using width 128.
  rescue_offset* {.intdefine: "somalier_q4_rescue_offset".} = 88

  ## Maximum physical span, in base pairs, from the first through last site of
  ## a rescue window. Wider physical windows are skipped even when they contain
  ## exactly `rescue_width` sites. Override with
  ## `-d:somalier_q4_rescue_max_span=10000000`.
  rescue_max_span* {.intdefine: "somalier_q4_rescue_max_span".} = 20_000_000

  ## Maximum allowed base-pair gap between adjacent sites in a rescue window.
  ## A window containing a larger gap is skipped so sparse marker deserts do
  ## not appear to provide continuous local compatibility. Override with
  ## `-d:somalier_q4_rescue_max_gap=1000000`.
  rescue_max_gap* {.intdefine: "somalier_q4_rescue_max_gap".} = 2_000_000

  ## Maximum count of jointly called opposite-homozygote sites (IBS0) allowed
  ## in a rescue window. Zero is the validated strict rule. Raising this value
  ## tolerates genotype errors but can sharply reduce specificity. Override
  ## with `-d:somalier_q4_rescue_max_ibs0=1`.
  rescue_max_ibs0* {.intdefine: "somalier_q4_rescue_max_ibs0".} = 0

  ## Minimum jointly called fraction of a rescue window, in permille. The
  ## derived count is `ceil(width * value / 1000)`; 500 at width 176 requires
  ## 88 jointly called sites. Override with
  ## `-d:somalier_q4_rescue_joint_permille=600`.
  rescue_joint_permille* {.intdefine: "somalier_q4_rescue_joint_permille".} = 500

  ## Minimum fraction of rescue-window sites that must be jointly callable and
  ## homozygous in each sample separately, in permille. At 250 and width 176,
  ## each sample needs 44 such homozygous sites. Override with
  ## `-d:somalier_q4_rescue_hom_permille=300`.
  rescue_hom_permille* {.intdefine: "somalier_q4_rescue_hom_permille".} = 250

  ## Minimum fraction of rescue-window sites with the same homozygous genotype
  ## in both samples, in permille. At 125 and width 176, 22 matching sites are
  ## required. Override with `-d:somalier_q4_rescue_match_permille=150`.
  rescue_match_permille* {.intdefine: "somalier_q4_rescue_match_permille".} = 125

  ## Typed projection seed consumed by the projection generator.
  projection_seed* = projection_seed_int.uint64
  ## Boolean form of `require_reciprocal_int` consumed by candidate generation.
  require_reciprocal* = require_reciprocal_int == 1
  ## Direct cosine threshold converted from thousandths to float32.
  direct_cosine_floor* = direct_cosine_milli.float32 / 1000'f32
  ## Rescue cosine threshold converted from thousandths to float32.
  rescue_cosine_floor* = rescue_cosine_milli.float32 / 1000'f32
  ## Quantization clipping limit converted from thousandths to float64.
  quantization_clip* = quantization_clip_milli.float64 / 1000.0
  ## Quantization multiplier converted from thousandths to float64.
  quantization_scale* = quantization_scale_milli.float64 / 1000.0
  ## Homozygous-reference boundary converted to float64 for production calling.
  hom_ref_cutoff* = hom_ref_cutoff_milli.float64 / 1000.0
  ## Homozygous-alternate boundary converted to float64 for production calling.
  hom_alt_cutoff* = hom_alt_cutoff_milli.float64 / 1000.0
  ## Minimum joint-call count, rounded up from width and permille.
  rescue_min_joint* =
    (rescue_width * rescue_joint_permille + 999) div 1000
  ## Minimum per-sample homozygous count, rounded up from width and permille.
  rescue_min_hom* =
    (rescue_width * rescue_hom_permille + 999) div 1000
  ## Minimum matching-homozygote count, rounded up from width and permille.
  rescue_min_match* =
    (rescue_width * rescue_match_permille + 999) div 1000

  ## Canonical text form of the mode, operating limits, and every
  ## method-defining compile-time value. Write this to logs. When persistent
  ## indexes are added, store the method-defining subset in metadata and
  ## require an exact match before loading an existing index. The diagnostic
  ## candidate-output path is omitted because it does not affect computation.
  q4_config_summary* =
    "mode=" & candidate_mode &
    ";exhaustive_max_samples=" & $exhaustive_max_samples &
    ";projection_batch_size=" & $projection_batch_size &
    ";threads=" & $query_threads &
    ";html_max_samples=" & $html_max_samples &
    ";dimensions=" & $projection_dimensions &
    ";seed=" & $projection_seed_int &
    ";m=" & $hnsw_connectivity &
    ";ef_construction=" & $hnsw_ef_construction &
    ";ef_search=" & $hnsw_ef_search &
    ";top_k=" & $neighbor_count &
    ";reciprocal=" & $require_reciprocal_int &
    ";direct_milli=" & $direct_cosine_milli &
    ";rescue_milli=" & $rescue_cosine_milli &
    ";clip_milli=" & $quantization_clip_milli &
    ";scale_milli=" & $quantization_scale_milli &
    ";hom_ref_cutoff_milli=" & $hom_ref_cutoff_milli &
    ";hom_alt_cutoff_milli=" & $hom_alt_cutoff_milli &
    ";rescue_width=" & $rescue_width &
    ";rescue_offset=" & $rescue_offset &
    ";rescue_max_span=" & $rescue_max_span &
    ";rescue_max_gap=" & $rescue_max_gap &
    ";rescue_max_ibs0=" & $rescue_max_ibs0 &
    ";rescue_joint_permille=" & $rescue_joint_permille &
    ";rescue_hom_permille=" & $rescue_hom_permille &
    ";rescue_match_permille=" & $rescue_match_permille

func use_q4_prefilter*(sample_count: int): bool {.inline.} =
  ## Resolve the compiled pair-generation policy for an expanded input count.
  case candidate_mode
  of "q4-hnsw": true
  of "exhaustive": false
  else: sample_count > exhaustive_max_samples

static:
  doAssert candidate_mode in ["auto", "exhaustive", "q4-hnsw"],
    "somalier_q4_candidate_mode must be auto, exhaustive, or q4-hnsw"
  doAssert exhaustive_max_samples >= 0,
    "somalier_q4_exhaustive_max_samples must be nonnegative"
  doAssert projection_batch_size > 0 and query_threads > 0,
    "Q4 projection batch size and thread count must be positive"
  doAssert html_max_samples >= 0,
    "somalier_q4_html_max_samples must be nonnegative"
  doAssert projection_dimensions > 0 and projection_dimensions mod 2 == 0,
    "usearchQ4Dimensions must be positive and even"
  doAssert projection_seed_int >= 0,
    "somalier_q4_projection_seed must be nonnegative"
  doAssert hnsw_connectivity > 0 and hnsw_ef_construction > 0 and
    hnsw_ef_search > 0 and neighbor_count > 0,
    "Q4 HNSW parameters must be positive"
  doAssert require_reciprocal_int in 0 .. 1,
    "somalier_q4_reciprocal must be 0 or 1"
  doAssert rescue_cosine_milli >= -1000 and direct_cosine_milli <= 1000 and
    rescue_cosine_milli < direct_cosine_milli,
    "Q4 cosine thresholds must be ordered and within cosine range"
  doAssert quantization_clip_milli > 0 and quantization_scale_milli > 0 and
    quantization_clip_milli * quantization_scale_milli <= 7_000_000,
    "Q4 clip times scale must fit signed levels -7 through +7"
  doAssert hom_ref_cutoff_milli >= 0 and hom_ref_cutoff_milli < 500 and
    hom_alt_cutoff_milli > 500 and hom_alt_cutoff_milli <= 1000 and
    hom_ref_cutoff_milli < hom_alt_cutoff_milli,
    "Q4 homozygous allele-balance cutoffs must straddle 0.5"
  doAssert rescue_width > 0 and rescue_offset >= 0 and
    rescue_offset < rescue_width,
    "Q4 rescue width and offset are inconsistent"
  doAssert rescue_max_span > 0 and rescue_max_gap > 0 and rescue_max_ibs0 >= 0,
    "Q4 rescue physical limits must be positive"
  doAssert rescue_joint_permille in 0 .. 1000 and
    rescue_hom_permille in 0 .. 1000 and rescue_match_permille in 0 .. 1000,
    "Q4 rescue support fractions must be between 0 and 1000 permille"
