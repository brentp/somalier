# Large-cohort relatedness with Q4/USearch

Somalier's default algorithm in N^2 in the number of samples.
Since it's so fast, this works for cohorts into the 10's of thousands.
Beyond that, it's too much.

When compiled with `-d:somalier_usearch`, `somalier` will use a sparse candidate
pre-filter search for large cohorts. This path projects each sample into a packed signed
4-bit vector, retrieves approximate neighbors with USearch HNSW, applies a
small positional rescue, and sends only admitted pairs to Somalier's exact
relatedness scorer.

The large-cohort implementation is excluded from ordinary builds. Enable it
with the `somalier_usearch` compile flag:

```bash
nimble build -d:somalier_usearch
# or directly with:
nim c -d:release -d:somalier_usearch src/somalier.nim
```

## Selecting exhaustive or sparse scoring


An ordinary build cannot execute the sparse path. If its compiled policy would
select Q4/HNSW, it exits with instructions to rebuild using
`-d:somalier_usearch`. To deliberately run exhaustive scoring above the
automatic boundary, compile with:

```bash
nimble build -d:somalier_q4_candidate_mode=exhaustive
```

That override permits quadratic memory and work and should be used only when
the cohort size is known to be manageable.

## Method

For every sample, Somalier calls genotypes using the same bitsets later read by
the exact pair scorer. Missing calls contribute zero to the projection and
neither match nor contradict during positional rescue. Consequently,
`-u/--unknown`, which converts missing calls to homozygous reference, is
rejected in Q4/HNSW mode.

The validated default method:

1. Standardizes called dosages using panel allele frequencies.
2. Projects them into 4,096 dimensions with the deterministic
   `splitmix64-rademacher-v1` projection and seed 1729.
3. Normalizes, clips, quantizes, and packs two signed 4-bit coordinates per
   byte. Each default sketch occupies 2,052 bytes including its squared norm.
4. Builds an in-memory HNSW index with `M=32` and `efConstruction=400`.
5. Queries 40 neighbors per sample with `efSearch=80` and keeps reciprocal
   neighbors.
6. Sends pairs with packed-Q4 cosine at least 0.20 directly to exact scoring.
7. For cosine from 0.11 through less than 0.20, admits pairs only when one of
   the two 176-site positional-window tilings passes the rescue rules.
8. Scores every admitted or pedigree-forced pair with Somalier's exact
   relatedness implementation.

The index and sketches are ephemeral. They are rebuilt on each invocation and
are not written to disk.

## Output behavior

Sparse mode writes the normal `samples.tsv`, `pairs.tsv`, and `groups.tsv`
files. `groups.tsv` remains an edge list of exact-scored pairs above its
existing relatedness threshold. Pedigree-forced pairs are added even when HNSW
does not retrieve them.

The exhaustive implementation randomly samples some unrelated pairs for
`pairs.tsv` and HTML. Sparse mode has no all-pairs stream from which to draw
that background, so it omits those random unrelated rows. It writes expected
pairs and candidates satisfying the existing output rule. HTML is skipped when
the sample count exceeds the compiled HTML limit because serializing every
sample would itself become large.

## Compile-time configuration

The authoritative documentation, validation checks, units, and interactions
for every setting are in
[`src/somalierpkg/usearch/q4_config.nim`](src/somalierpkg/usearch/q4_config.nim).
Numeric fractional settings use integer thousandths or permille because Nim
compile-time defines are integers.

| Compile define | Default | Purpose |
|---|---:|---|
| `somalier_q4_candidate_mode` | `auto` | Select `auto`, `exhaustive`, or `q4-hnsw`. Q4 execution still requires `-d:somalier_usearch`. |
| `somalier_q4_exhaustive_max_samples` | `5000` | Largest exhaustively scored cohort under `auto`; Q4 starts above this count. |
| `somalier_q4_projection_batch_size` | `256` | Samples per projection batch; trades working memory for throughput without changing candidates. |
| `somalier_q4_threads` | `1` | Concurrent HNSW queries and reserved search slots. |
| `somalier_q4_html_max_samples` | `10000` | Largest sparse cohort receiving HTML; zero disables sparse HTML. |
| `somalier_q4_candidate_output` | empty | Optional TSV path for every admitted pair, Q4 cosine, and direct/rescue reason. The file is overwritten each run. |
| `usearchQ4Dimensions` | `4096` | Projection coordinates and packed-sketch size; must be positive and even. |
| `somalier_q4_projection_seed` | `1729` | Seed for the deterministic Rademacher projection. |
| `somalier_q4_m` | `32` | HNSW graph connectivity; affects index memory, construction, and retrieval. |
| `somalier_q4_ef_construction` | `400` | HNSW candidate breadth while constructing the graph. |
| `somalier_q4_ef_search` | `80` | HNSW expansion during neighbor queries; normally at least top K. |
| `somalier_q4_top_k` | `40` | Non-self neighbors requested for each sample. |
| `somalier_q4_reciprocal` | `1` | `1` requires mutual retrieval; `0` uses the union of directed neighbors. |
| `somalier_q4_direct_cosine_milli` | `200` | Q4 cosine threshold for direct exact scoring; `200` means 0.200. |
| `somalier_q4_rescue_cosine_milli` | `110` | Lower cosine boundary for positional rescue; pairs below it are discarded. |
| `somalier_q4_clip_milli` | `3500` | Absolute normalized-coordinate clipping limit; `3500` means 3.5. |
| `somalier_q4_scale_milli` | `2000` | Multiplier before rounding to signed Q4; `2000` means 2.0. |
| `somalier_q4_hom_ref_cutoff_milli` | `10` | Strict hom-ref allele-balance upper boundary; `10` means AB below 0.01. |
| `somalier_q4_hom_alt_cutoff_milli` | `990` | Strict hom-alt allele-balance lower boundary; `990` means AB above 0.99. |
| `somalier_q4_rescue_width` | `176` | Ordered panel sites in each positional-rescue window. |
| `somalier_q4_rescue_offset` | `88` | Starting offset for the second rescue-window tiling. |
| `somalier_q4_rescue_max_span` | `20000000` | Maximum physical span of a rescue window in base pairs. |
| `somalier_q4_rescue_max_gap` | `2000000` | Maximum adjacent-site gap inside a rescue window in base pairs. |
| `somalier_q4_rescue_max_ibs0` | `0` | Maximum opposite-homozygote sites allowed in a rescue window. |
| `somalier_q4_rescue_joint_permille` | `500` | Minimum jointly called fraction of a rescue window; 500 means 50%. |
| `somalier_q4_rescue_hom_permille` | `250` | Minimum per-sample jointly callable homozygous fraction; 250 means 25%. |
| `somalier_q4_rescue_match_permille` | `125` | Minimum matching-homozygote fraction; 125 means 12.5%. |

Multiple defines can be supplied to one build:

```bash
nimble build -d:somalier_usearch \
  -d:usearchQ4Dimensions=4096 \
  -d:somalier_q4_threads=8 \
  -d:somalier_q4_direct_cosine_milli=220 \
  -d:somalier_q4_rescue_cosine_milli=115
```

Every Q4 run logs the complete effective configuration.

## Parameters useful for research

The most useful parameters to vary are:

1. **`somalier_q4_direct_cosine_milli`:** pairs at or above this cosine go
   directly to exact scoring. Lowering it improves recall by bypassing rescue
   for more pairs, but increases exact-scoring work. Raising it moves pairs
   into the rescue band, where they can be rejected.
2. **`somalier_q4_rescue_cosine_milli`:** pairs below this cosine are discarded;
   pairs from here to the direct threshold must pass positional rescue.
   Lowering it tests more weak neighbors and can recover relatives at the cost
   of more rescue work and candidates.
3. **`somalier_q4_top_k`:** controls how many neighbors each sample retrieves
   before reciprocity and cosine filtering. Larger K can recover relatives
   hidden in dense neighborhoods, but increases query storage and the maximum
   reciprocal candidate count (`N*K/2`).
4. **`somalier_q4_reciprocal`:** `1` requires both samples to retrieve each
   other and is the main specificity control. `0` accepts the directed-neighbor
   union, usually improving recall while substantially increasing candidates.
5. **`usearchQ4Dimensions`:** controls projection capacity and packed-vector
   size. More dimensions may preserve relatedness signal better, while index
   memory, projection work, and vector-distance work grow approximately
   linearly.
6. **`somalier_q4_projection_seed`:** changes the deterministic random
   projection without changing its size. Sweep several fixed seeds and report
   worst-seed performance; selecting the best seed on one cohort risks
   overfitting.
7. **`somalier_q4_ef_search`:** controls how broadly HNSW explores the graph at
   query time. Increasing it can recover true top-K neighbors missed by the
   approximate search, but does not change K and costs query time.
8. **`somalier_q4_m` and `somalier_q4_ef_construction`:** M controls graph
   connectivity and persistent index memory; `ef_construction` controls search
   breadth while building that graph. Increasing either can improve graph
   quality, with M costing more memory and both costing build time.
9. **`somalier_q4_rescue_width` and `somalier_q4_rescue_offset`:** width is the
   number of ordered panel sites per rescue window; offset starts the second
   tiling. Wider windows are generally more specific but more sensitive to
   missing calls and errors. When changing width, normally keep the offset near
   half the width so the second tiling covers first-tiling boundaries.
10. **`somalier_q4_rescue_max_span` and `somalier_q4_rescue_max_gap`:** reject
    windows that cover too much physical distance or contain a large marker
    gap. Lower values demand stronger local continuity, but can remove most
    windows on sparse panels.
11. **`somalier_q4_rescue_max_ibs0`:** limits opposite-homozygote sites within
    a passing window. Raising it tolerates genotype errors, but weakens one of
    rescue's strongest unrelated-pair filters and can sharply reduce
    specificity.
12. **`somalier_q4_rescue_joint_permille`,
    `somalier_q4_rescue_hom_permille`, and
    `somalier_q4_rescue_match_permille`:** require joint calls, homozygous calls
    in each sample, and matching homozygotes, respectively. Counts are rounded
    up from width times permille, so width and these fractions must be evaluated
    together. Higher values increase evidence per window but penalize missing
    or noisy samples.
13. **`somalier_q4_clip_milli` and `somalier_q4_scale_milli`:** jointly map
    normalized projections into signed values from -7 through +7. More scale
    gives finer resolution near zero but can saturate tails; clipping limits
    those tails. Their product may not exceed 7 after unit conversion, and any
    change alters cosine geometry.
14. **`somalier_q4_hom_ref_cutoff_milli` and
    `somalier_q4_hom_alt_cutoff_milli`:** set the strict allele-balance limits
    for homozygous calls. They affect projection, rescue, and exact pair
    statistics, so treat them as caller changes and validate more broadly than
    an HNSW-only adjustment.

`somalier_q4_candidate_output` is especially useful during these experiments:
it preserves candidate IDs, exact packed-Q4 cosine values, and admission
reasons. Compare candidate recall separately from recall after the exact
relatedness threshold.

The remaining settings are operational controls.
`somalier_q4_projection_batch_size` trades temporary projection memory for
throughput, and `somalier_q4_threads` controls query concurrency; neither
should change candidates. `somalier_q4_html_max_samples` changes only whether
HTML is written. `somalier_q4_exhaustive_max_samples` and
`somalier_q4_candidate_mode` select the execution path and are useful for
integration and exhaustive-parity testing, not for improving retrieval.
