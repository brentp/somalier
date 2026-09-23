import unittest

import somalierpkg/usearch/q4_config

suite "Q4 compile-time defaults":
  test "match the validated profile":
    check candidate_mode == "auto"
    check exhaustive_max_samples == 5000
    check projection_batch_size == 256
    check query_threads == 1
    check html_max_samples == 10_000
    check candidate_output_path == ""
    check not use_q4_prefilter(5000)
    check use_q4_prefilter(5001)
    check projection_dimensions == 4096
    check projection_seed == 1729'u64
    check hnsw_connectivity == 32
    check hnsw_ef_construction == 400
    check hnsw_ef_search == 80
    check neighbor_count == 40
    check require_reciprocal
    check direct_cosine_floor == 0.20'f32
    check rescue_cosine_floor == 0.11'f32
    check quantization_clip == 3.5
    check quantization_scale == 2.0
    check hom_ref_cutoff == 0.01
    check hom_alt_cutoff == 0.99
    check rescue_width == 176
    check rescue_offset == 88
    check rescue_min_joint == 88
    check rescue_min_hom == 44
    check rescue_min_match == 22
    check q4_config_summary.len > 0
