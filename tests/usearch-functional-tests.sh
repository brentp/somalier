#!/usr/bin/env bash

set -euo pipefail

workdir=$(mktemp -d "${TMPDIR:-/tmp}/somalier-usearch-test.XXXXXX")
trap 'rm -rf "$workdir"' EXIT

baseline=./src/somalier
usearch_binary="$workdir/somalier-usearch"
candidate_manifest="$workdir/q4-candidates.tsv"

if [[ ! -x "$baseline" ]]; then
  echo "baseline binary is missing; run tests/functional-tests.sh first" >&2
  exit 1
fi

nim c \
  -d:debug \
  -d:useSysAssert \
  -d:useGcAssert \
  -d:somalier_usearch \
  -d:somalier_q4_candidate_mode=q4-hnsw \
  "-d:somalier_q4_candidate_output=$candidate_manifest" \
  --lineDir:on \
  --debuginfo \
  --boundChecks:on \
  -x:on \
  "--nimcache:$workdir/nimcache" \
  "--out:$usearch_binary" \
  src/somalier

"$usearch_binary" extract \
  --sample-prefix Q4-A- \
  --sites tests/test_sites.vcf \
  --fasta tests/test.fa \
  --out-dir "$workdir/a" \
  tests/gt_only.vcf.gz

"$usearch_binary" extract \
  --sample-prefix Q4-B- \
  --sites tests/test_sites.vcf \
  --fasta tests/test.fa \
  --out-dir "$workdir/b" \
  tests/gt_only.vcf.gz

mapfile -t extracted < <(find "$workdir/a" "$workdir/b" -name '*.somalier' -type f | sort)
if [[ ${#extracted[@]} -ne 2 ]]; then
  echo "expected two extracted samples, found ${#extracted[@]}" >&2
  exit 1
fi

"$baseline" relate \
  --sites tests/test_sites.vcf \
  --output-prefix "$workdir/exhaustive" \
  "${extracted[@]}"

if ! "$usearch_binary" relate \
  --sites tests/test_sites.vcf \
  --output-prefix "$workdir/q4" \
  "${extracted[@]}" \
  2>"$workdir/q4.stderr"; then
  cat "$workdir/q4.stderr" >&2
  exit 1
fi
cat "$workdir/q4.stderr" >&2

grep -Fq "Q4 candidate mode is experimental" "$workdir/q4.stderr"
grep -Fq "mode=q4-hnsw" "$workdir/q4.stderr"
grep -Fq "Q4 candidates reciprocal=" "$workdir/q4.stderr"

test "$(wc -l < "$candidate_manifest")" -eq 2
grep -Eq $'Q4-A-test_sample\tQ4-B-test_sample\t.*\tdirect$' "$candidate_manifest"

test "$(wc -l < "$workdir/q4.samples.tsv")" -eq 3
test "$(wc -l < "$workdir/q4.pairs.tsv")" -eq 2
diff -u "$workdir/exhaustive.pairs.tsv" "$workdir/q4.pairs.tsv"

echo "USearch functional tests passed"
