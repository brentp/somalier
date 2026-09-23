## Q4 HNSW candidate generation and sparse relate output.
##
## This file is included by relate.nim so it can reuse the established private
## reporting and pedigree helpers without exporting or moving the exhaustive
## implementation. All USearch-specific code remains in this directory.

import arraymancer
import usearch_q4
import std/atomics

static:
  doAssert compileOption("threads"),
    "somalier_usearch requires Nim thread support"

const
  q4_projection_version = "splitmix64-rademacher-v1"

type
  q4_rescue_window = object
    start, stop: int

  q4_load_result = object
    final: relation_matrices
    charr_stats: seq[CharrStats]
    candidates: seq[uint64]
    direct_count: int
    rescued_count: int
    reciprocal_count: int

  q4_query_job = object
    index: ptr Q4Index
    directed: ptr UncheckedArray[uint32]
    counts: ptr UncheckedArray[int]
    next_sample: Atomic[int]
    sample_count: int
    wanted: int

when compileOption("threads"):
  proc q4_query_worker(job: ptr q4_query_job) {.thread.} =
    while true:
      let sample = job.next_sample.fetchAdd(1)
      if sample >= job.sample_count: break
      let found = job.index[].searchById(sample.uint32, job.wanted)
      job.counts[sample] = found.len
      for i, neighbor in found:
        job.directed[sample * job.wanted + i] = neighbor.id

proc q4_read_extracted(path: string, cnt: var counts) =
  var stream = newFileStream(path, fmRead)
  if stream == nil:
    raise newException(IOError, "could not open extracted file: " & path)
  defer: stream.close()
  var name_length: uint8
  var version: uint8
  discard stream.readData(version.addr, sizeof(version))
  if version != formatVersion:
    raise newException(ValueError, &"expected format version {formatVersion}, got {version} in {path}")
  discard stream.readData(name_length.addr, sizeof(name_length))
  cnt.sample_name = newString(name_length.int)
  if name_length > 0:
    discard stream.readData(cnt.sample_name[0].addr, name_length.int)
  var autosomal_count, x_count, y_count: uint16
  discard stream.readData(autosomal_count.addr, sizeof(autosomal_count))
  discard stream.readData(x_count.addr, sizeof(x_count))
  discard stream.readData(y_count.addr, sizeof(y_count))
  cnt.sites = newSeq[allele_count](autosomal_count.int)
  cnt.x_sites = newSeq[allele_count](x_count.int)
  cnt.y_sites = newSeq[allele_count](y_count.int)
  for values in [addr cnt.sites, addr cnt.x_sites, addr cnt.y_sites]:
    if values[].len > 0:
      let wanted = values[].len * sizeof(allele_count)
      if stream.readData(values[][0].addr, wanted) != wanted:
        raise newException(IOError, "truncated extracted file: " & path)

proc next_split_mix64(state: var uint64): uint64 {.inline.} =
  state += 0x9e3779b97f4a7c15'u64
  var z = state
  z = (z xor (z shr 30)) * 0xbf58476d1ce4e5b9'u64
  z = (z xor (z shr 27)) * 0x94d049bb133111eb'u64
  z xor (z shr 31)

proc q4_alts(balance, min_ab: float64): int8 {.inline.} =
  if balance < 0: return -1
  if balance < hom_ref_cutoff: return 0
  if balance > hom_alt_cutoff: return 2
  if balance >= min_ab and balance <= 1.0 - min_ab: return 1
  -1

proc q4_autosomal_call(c: allele_count, min_ab: float64,
    min_depth: int): tuple[balance: float64, genotype: int8] {.inline.} =
  result.balance = c.ab(min_depth)
  result.genotype = result.balance.q4_alts(min_ab)

proc q4_x_call(c: allele_count, min_ab: float64,
    min_depth: int): int8 {.inline.} =
  let depth = c.nref + c.nalt
  let total = depth + c.nother
  if total > 0 and c.nother.float64 / total.float64 > 0.04: return -1
  if depth.int < min_depth: return -1
  let balance = if depth == 0: -1.0 else: c.nalt.float64 / depth.float64
  balance.q4_alts(min_ab)

proc set_arena_bit(arena: var seq[uint64], sample, words, genotype,
    site: int) {.inline.} =
  let offset = sample * 3 * words + genotype * words + (site shr 6)
  arena[offset] = arena[offset] or (1'u64 shl (site and 63))

proc arena_call(final: relation_matrices, sample, site: int): int8 {.inline.} =
  let word = site shr 6
  let bit = 1'u64 shl (site and 63)
  let base = sample * 3 * final.autosomal_words
  if (final.autosomal_arena[base + word] and bit) != 0: return 0
  if (final.autosomal_arena[base + final.autosomal_words + word] and bit) != 0: return 1
  if (final.autosomal_arena[base + 2 * final.autosomal_words + word] and bit) != 0: return 2
  -1

proc q4_projection_matrix(site_count: int): Tensor[float32] =
  result = newTensorUninit[float32](site_count, projection_dimensions)
  var state = projection_seed
  for site in 0 ..< site_count:
    for coordinate in 0 ..< projection_dimensions:
      result[site, coordinate] =
        if (next_split_mix64(state) and 1'u64) == 0: -1'f32 else: 1'f32

proc q4_standardized_batch(final: relation_matrices,
    frequencies: openArray[float32], first, stop: int): Tensor[float32] =
  result = newTensor[float32](stop - first, frequencies.len)
  for sample in first ..< stop:
    for site, p in frequencies:
      let genotype = final.arena_call(sample, site)
      if genotype >= 0:
        let variance = max(2'f32 * p * (1'f32 - p), 1e-6'f32)
        result[sample - first, site] =
          (genotype.float32 - 2'f32 * p) / sqrt(variance)

proc q4_quantize(projected: Tensor[float32], sample: int): Q4Record =
  var norm2 = 0'f64
  for coordinate in 0 ..< projection_dimensions:
    let value = projected[sample, coordinate].float64
    norm2 += value * value
  if norm2 == 0:
    raise newException(ValueError, "Q4 projection has zero norm")
  let normalization = sqrt(projection_dimensions.float64 / norm2)
  var values: array[projection_dimensions, int8]
  for coordinate in 0 ..< projection_dimensions:
    let normalized = projected[sample, coordinate].float64 * normalization
    let clipped = clamp(normalized, -quantization_clip, quantization_clip)
    values[coordinate] = round(clipped * quantization_scale).int8
  packQ4(values)

proc q4_rescue_windows(sites: openArray[Site]): seq[q4_rescue_window] =
  var chrom_start = 0
  while chrom_start < sites.len:
    var chrom_stop = chrom_start + 1
    while chrom_stop < sites.len and
        sites[chrom_stop].chrom == sites[chrom_start].chrom:
      chrom_stop.inc
    for offset in [0, rescue_offset]:
      var start = chrom_start + offset
      while start + rescue_width <= chrom_stop:
        let stop = start + rescue_width
        var max_gap = 0
        for i in start + 1 ..< stop:
          max_gap = max(max_gap, sites[i].position - sites[i - 1].position)
        let span = sites[stop - 1].position - sites[start].position
        if span <= rescue_max_span and max_gap <= rescue_max_gap:
          result.add(q4_rescue_window(start: start, stop: stop))
        start += rescue_width
    chrom_start = chrom_stop

proc q4_passes_rescue(final: relation_matrices, first, second: int,
    windows: openArray[q4_rescue_window]): bool =
  for window in windows:
    var joint, first_hom, second_hom, matching_hom, opposite: int
    for site in window.start ..< window.stop:
      let a = final.arena_call(first, site)
      let b = final.arena_call(second, site)
      if a < 0 or b < 0: continue
      joint.inc
      if a != 1: first_hom.inc
      if b != 1: second_hom.inc
      if a != 1 and b != 1:
        if a == b: matching_hom.inc
        else: opposite.inc
    if opposite <= rescue_max_ibs0 and joint >= rescue_min_joint and
        first_hom >= rescue_min_hom and second_hom >= rescue_min_hom and
        matching_hom >= rescue_min_match:
      return true

proc q4_contains(neighbors: openArray[uint32], start, count: int,
    id: uint32): bool {.inline.} =
  for i in start ..< start + count:
    if neighbors[i] == id: return true

proc pack_pair(first, second: int): uint64 {.inline.} =
  (first.uint64 shl 32) or second.uint32.uint64

proc unpack_pair(value: uint64): tuple[first, second: int] {.inline.} =
  ((value shr 32).int, (value and 0xffff_ffff'u64).int)

proc q4_load_and_index(paths: seq[string], sites_path: string,
    min_ab: float64, min_depth, batch_size, threads: int,
    charr_hom_rate, charr_hom_alpha: float64,
    candidate_output_path = ""): q4_load_result =
  let load_started = epochTime()
  let panel = readSitesWithAF(sites_path)
  let all_sites = readSites(sites_path)
  var flips, x_flips: seq[bool]
  var expected_y_sites = 0
  for site in all_sites:
    case site.chrom
    of "X", "chrX", "NC_000023.10", "NC_000023.11": x_flips.add(site.flip)
    of "Y", "chrY", "NC_000024.9", "NC_000024.10": expected_y_sites.inc
    else: flips.add(site.flip)

  var ordered_paths = paths.sorted
  let sample_count = ordered_paths.len
  result.final.samples = newSeq[string](sample_count)
  result.final.stats = newSeq[Stat4](sample_count)
  result.charr_stats = newSeq[CharrStats](sample_count)
  for genotype in 0 ..< result.final.gt_counts.len:
    result.final.gt_counts[genotype] = newSeq[uint16](sample_count)

  var cnt: counts
  var seen = initHashSet[string]()
  for sample, path in ordered_paths:
    q4_read_extracted(path, cnt)
    if cnt.sites.len != panel.sites.len or cnt.x_sites.len != x_flips.len or
        cnt.y_sites.len != expected_y_sites:
      raise newException(ValueError, &"[somalier] extracted input {path} has " &
        &"{cnt.sites.len}/{cnt.x_sites.len}/{cnt.y_sites.len} autosomal/X/Y sites; " &
        &"expected {panel.sites.len}/{x_flips.len}/{expected_y_sites} from --sites")
    if cnt.sample_name in seen:
      raise newException(ValueError, "duplicate sample name: " & cnt.sample_name)
    seen.incl(cnt.sample_name)
    result.final.samples[sample] = cnt.sample_name
    if sample == 0:
      result.final.autosomal_words = (cnt.sites.len + 63) div 64
      result.final.x_words = (cnt.x_sites.len + 63) div 64
      result.final.autosomal_arena = newSeq[uint64](sample_count * 3 *
        result.final.autosomal_words)
      result.final.x_arena = newSeq[uint64](sample_count * 3 *
        result.final.x_words)

    var stat: Stat4
    for site, c in cnt.sites:
      let called = q4_autosomal_call(c, min_ab, min_depth)
      if flips[site] and called.genotype >= 0:
        stat.hom_ref_adjustment += called.genotype.int32 - 1
      stat.dp.push(int(c.nref + c.nalt))
      if c.nref > 0 or c.nalt > 0 or c.nother > 0:
        stat.un.push(c.nother.float64 / (c.nref + c.nalt + c.nother).float64)
      if c.nref.float64 > min_depth.float64 / 2 or
          c.nalt.float64 > min_depth.float64 / 2:
        stat.ab.push(called.balance)
      if called.balance != -1: stat.gtdp.push(int(c.nref + c.nalt))
      if called.balance > 0.02 and called.balance < 0.98 and
          (called.balance < 0.1 or called.balance > 0.9):
        result.final.gt_counts[4][sample].inc
      if called.genotype < 0:
        result.final.gt_counts[3][sample].inc
      else:
        result.final.gt_counts[called.genotype][sample].inc
        result.final.autosomal_arena.set_arena_bit(sample,
          result.final.autosomal_words, called.genotype, site)

    for site, c in cnt.x_sites:
      let genotype = q4_x_call(c, min_ab, min_depth)
      if genotype < 0: continue
      stat.x_dp.push((c.nref + c.nalt).float)
      if x_flips[site]:
        stat.x_hom_ref_adjustment += genotype.int32 - 1
      case genotype
      of 0:
        stat.x_hom_ref.inc
      of 1:
        stat.x_het.inc
      of 2:
        stat.x_hom_alt.inc
      else: discard
      result.final.x_arena.set_arena_bit(sample, result.final.x_words,
        genotype, site)
    for c in cnt.y_sites:
      if q4_x_call(c, min_ab, min_depth) >= 0:
        stat.y_dp.push((c.nref + c.nalt).float)
    result.final.stats[sample] = stat
    result.charr_stats[sample] = estimate_charr(cnt.sites, panel.pop_afs,
      min_depth, charr_hom_rate, charr_hom_alpha)
  stderr.write_line &"[somalier] Q4 read and called {sample_count} samples in {epochTime() - load_started:.2f}s"

  let projection_started = epochTime()
  var projection = q4_projection_matrix(panel.sites.len)
  let requested_query_threads = min(threads, max(1, sample_count))
  var index = createQ4Index(sample_count, connectivity = hnsw_connectivity,
    expansionAdd = hnsw_ef_construction, expansionSearch = hnsw_ef_search,
    searchThreads = requested_query_threads)
  # USearch fails rather than waits when all search slots are occupied. Read
  # back the effective capacity so a wrapper/backend adjustment cannot cause
  # Somalier to launch more concurrent searches than the index can serve.
  let search_slots = index.searchThreads
  if search_slots < 1:
    raise newException(ValueError, "USearch index has no query slots")
  let query_workers = min(requested_query_threads, search_slots)
  var first = 0
  while first < sample_count:
    let stop = min(sample_count, first + batch_size)
    var values = q4_standardized_batch(result.final, panel.pop_afs, first, stop)
    var projected = values * projection
    values = Tensor[float32]()
    for sample in first ..< stop:
      discard index.add(q4_quantize(projected, sample - first))
    projected = Tensor[float32]()
    first = stop
  projection = Tensor[float32]()
  stderr.write_line &"[somalier] Q4 projected and indexed in {epochTime() - projection_started:.2f}s (batch={batch_size})"

  if sample_count < 2: return
  let query_started = epochTime()
  let wanted = min(neighbor_count, sample_count - 1)
  var directed = newSeq[uint32](sample_count * wanted)
  var counts = newSeq[int](sample_count)
  when compileOption("threads"):
    var job = q4_query_job(
      index: addr index,
      directed: cast[ptr UncheckedArray[uint32]](addr directed[0]),
      counts: cast[ptr UncheckedArray[int]](addr counts[0]),
      sample_count: sample_count,
      wanted: wanted)
    job.next_sample.store(0)
    var workers = newSeq[Thread[ptr q4_query_job]](query_workers)
    for worker in workers.mitems:
      createThread(worker, q4_query_worker, addr job)
    joinThreads(workers)
  else:
    for sample in 0 ..< sample_count:
      let found = index.searchById(sample.uint32, wanted)
      counts[sample] = found.len
      for i, neighbor in found:
        directed[sample * wanted + i] = neighbor.id

  let windows = q4_rescue_windows(panel.sites)
  var candidate_output: File
  if candidate_output_path.len > 0:
    if not open(candidate_output, candidate_output_path, fmWrite):
      raise newException(IOError, "could not open Q4 candidate output: " &
        candidate_output_path)
    candidate_output.write_line "sample_a\tsample_b\tq4_cosine\tadmission"
  for first in 0 ..< sample_count:
    for i in 0 ..< counts[first]:
      let second = directed[first * wanted + i].int
      let reciprocal = q4_contains(directed, second * wanted, counts[second],
        first.uint32)
      if require_reciprocal:
        if second <= first or not reciprocal: continue
      elif reciprocal and second < first:
        continue
      let low = min(first, second)
      let high = max(first, second)
      result.reciprocal_count.inc
      let cosine = index.cosineById(low.uint32, high.uint32)
      if cosine >= direct_cosine_floor:
        result.direct_count.inc
        result.candidates.add(pack_pair(low, high))
        if candidate_output_path.len > 0:
          candidate_output.write_line &"{result.final.samples[low]}\t{result.final.samples[high]}\t{cosine:.7f}\tdirect"
      elif cosine >= rescue_cosine_floor and
          q4_passes_rescue(result.final, low, high, windows):
        result.rescued_count.inc
        result.candidates.add(pack_pair(low, high))
        if candidate_output_path.len > 0:
          candidate_output.write_line &"{result.final.samples[low]}\t{result.final.samples[high]}\t{cosine:.7f}\trescue"
  if candidate_output_path.len > 0: candidate_output.close()
  if candidate_output_path.len > 0:
    stderr.write_line "[somalier] wrote Q4 candidate manifest to: " &
      candidate_output_path
  result.candidates.sort
  stderr.write_line &"[somalier] Q4 queried with {query_workers} thread(s) and filtered candidates in {epochTime() - query_started:.2f}s"

proc add_forced_pairs(candidates: var seq[uint64], groups: openArray[pair],
    sample_names: openArray[string]) =
  var by_name = initTable[string, int]()
  for i, name in sample_names: by_name[name] = i
  for group in groups:
    if group.a notin by_name or group.b notin by_name: continue
    let first = min(by_name[group.a], by_name[group.b])
    let second = max(by_name[group.a], by_name[group.b])
    if first != second: candidates.add(pack_pair(first, second))
  candidates.sort
  if candidates.len > 1:
    var write_at = 1
    for read_at in 1 ..< candidates.len:
      if candidates[read_at] != candidates[write_at - 1]:
        candidates[write_at] = candidates[read_at]
        write_at.inc
    candidates.setLen(write_at)

proc run_q4_relate(paths: seq[string], sites_path, groups_path,
    ped_path, output_prefix: string, sample_prefix: seq[string],
    min_depth: int, min_ab, charr_hom_rate, charr_hom_alpha: float64,
    infer: bool) =
  stderr.write_line "[somalier] Q4 candidate mode is experimental"
  stderr.write_line "[somalier] Q4 compile_config " & q4_config_summary
  stderr.write_line "[somalier] Q4 projection=" & q4_projection_version &
    " metric=" & metricImplementation()
  stderr.write_line "[somalier] Q4 validates panel site counts; the current .somalier format does not store a panel checksum"
  stderr.write_line &"[somalier] Q4 caller min_depth={min_depth} min_ab={min_ab} " &
    &"hom_ref_ab<{hom_ref_cutoff} hom_alt_ab>{hom_alt_cutoff}"
  if min_depth != 7 or abs(min_ab - 0.3) > 1e-9 or
      abs(hom_ref_cutoff - 0.01) > 1e-12 or
      abs(hom_alt_cutoff - 0.99) > 1e-12:
    stderr.write_line "[somalier] WARNING: Q4 retrieval was validated with --min-depth 7, --min-ab 0.3, hom-ref AB < 0.01, and hom-alt AB > 0.99"

  var t0 = cpuTime()
  var loaded = q4_load_and_index(paths, sites_path, min_ab, min_depth,
    projection_batch_size, query_threads, charr_hom_rate, charr_hom_alpha,
    candidate_output_path)
  stderr.write_line &"[somalier] Q4 loaded, projected, and indexed {loaded.final.samples.len} samples in {cpuTime() - t0:.2f}s"

  var samples: seq[Sample]
  var groups: seq[pair]
  if ped_path != "": samples = parse_ped(ped_path)
  groups.add_ped_samples(samples, loaded.final.samples)
  groups.add_prefixed_samples(loaded.final.samples, sample_prefix)
  groups.add(readGroups(groups_path, groups))
  groups.sort(cmp_pair)
  loaded.candidates.add_forced_pairs(groups, loaded.final.samples)

  var fh_tsv, fh_samples: File
  if not open(fh_tsv, output_prefix & "pairs.tsv", fmWrite):
    quit "couldn't open output file"
  if not open(fh_samples, output_prefix & "samples.tsv", fmWrite):
    quit "couldn't open output file"
  fh_tsv.write_line '#', header.replace("$", "")

  let write_html = html_max_samples > 0 and
    loaded.final.samples.len <= html_max_samples
  var fh_html: File
  var tmpls: seq[string]
  var rels: seq[relations]
  if write_html:
    if not open(fh_html, output_prefix & "html", fmWrite):
      quit "couldn't open html output file"
    tmpls = tmpl_html.split("<INPUT_JSON>")
    fh_html.write(tmpls[0].replace("<SAMPLE_JSON>", toj(loaded.final.samples,
      loaded.final.stats, loaded.final.gt_counts, loaded.charr_stats,
      samples.to_sex_lookup)))
  else:
    stderr.write_line &"[somalier] skipping HTML for {loaded.final.samples.len} samples (compiled Q4 HTML limit={html_max_samples})"

  let check_hom_alts = getEnv("SOMALIER_CHECK_HOM_ALTS") != ""
  var grouped: seq[pair]
  var parent_child_pair = newTable[string, seq[string]]()
  var sib_pairs = newTable[string, seq[string]]()
  var rel_gt_0p2 = newTable[string, seq[string]]()
  var written_pairs = 0
  for packed in loaded.candidates:
    let ids = unpack_pair(packed)
    var rel = loaded.final.relatedness(ids.first, ids.second)
    rel.sample_a = loaded.final.samples[ids.first]
    rel.sample_b = loaded.final.samples[ids.second]
    if rel.rel > 0.125:
      grouped.add((rel.sample_a, rel.sample_b, rel.rel))
    var group_index = groups.binarySearch((rel.sample_a, rel.sample_b, -1.0),
      cmp_pair)
    if group_index == -1:
      group_index = groups.binarySearch((rel.sample_b, rel.sample_a, -1.0),
        cmp_pair)
    let expected = if group_index == -1: -1.0 else: groups[group_index].rel
    let rr = rel.rel
    let ibs_ratio = rel.ibs0.float / rel.ibs2.float
    if rr > 0.4 and rr < 0.6 and ibs_ratio < 0.005:
      parent_child_pair.mgetOrPut(rel.sample_a, @[]).add(rel.sample_b)
      parent_child_pair.mgetOrPut(rel.sample_b, @[]).add(rel.sample_a)
    elif rr > 0.38 and rr < 0.62 and ibs_ratio > 0.015 and ibs_ratio < 0.052:
      sib_pairs.mgetOrPut(rel.sample_a, @[]).add(rel.sample_b)
      sib_pairs.mgetOrPut(rel.sample_b, @[]).add(rel.sample_a)
    elif rr > 0.96 and ibs_ratio < 0.005:
      var specified_as_sibs = false
      for group in groups:
        if ((group.a == rel.sample_a and group.b == rel.sample_b) or
            (group.b == rel.sample_a and group.a == rel.sample_b)) and
            group.rel > 0.4:
          specified_as_sibs = true
          break
      if specified_as_sibs and infer:
        sib_pairs.mgetOrPut(rel.sample_a, @[]).add(rel.sample_b)
        sib_pairs.mgetOrPut(rel.sample_b, @[]).add(rel.sample_a)
    elif rr > 0.2 and loaded.final.gt_counts.high_quality(ids.first,
        check_hom_alts) and loaded.final.gt_counts.high_quality(ids.second,
        check_hom_alts):
      rel_gt_0p2.mgetOrPut(rel.sample_a, @[]).add(rel.sample_b)
      rel_gt_0p2.mgetOrPut(rel.sample_b, @[]).add(rel.sample_a)

    let interesting = expected != -1 or rr > 0.05
    if interesting:
      fh_tsv.write_line rel.tsv(expected)
      written_pairs.inc
      if write_html: rels.add(rel, max(0, expected))

  if write_html:
    fh_html.write(%* rels)
    fh_html.write(tmpls[1])
    fh_html.close()
    stderr.write_line "[somalier] wrote Q4 candidate HTML to: " & output_prefix & "html"

  var looker = loaded.final.look(samples, loaded.final.stats,
    parent_child_pair, sib_pairs, rel_gt_0p2)
  fh_samples.write_ped(loaded.final, loaded.final.stats,
    loaded.final.gt_counts, loaded.charr_stats, looker, infer)
  fh_tsv.close()
  grouped.write(output_prefix)
  stderr.write_line &"[somalier] Q4 candidates reciprocal={loaded.reciprocal_count} direct={loaded.direct_count} rescued={loaded.rescued_count} scored={loaded.candidates.len} written={written_pairs} possible={loaded.final.samples.len.int64 * (loaded.final.samples.len - 1).int64 div 2}"
  stderr.write_line "[somalier] wrote groups to: " & output_prefix & "groups.tsv"
  stderr.write_line "[somalier] wrote samples to: " & output_prefix & "samples.tsv"
  stderr.write_line "[somalier] wrote pair-wise relatedness metrics to: " &
    output_prefix & "pairs.tsv"
