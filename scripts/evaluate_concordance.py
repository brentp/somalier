#!/usr/bin/env python3
import argparse
import csv
import statistics
from pathlib import Path


def auc(scores):
    scores = sorted(scores)
    n1 = sum(label for _, label in scores)
    n0 = len(scores) - n1
    if n1 == 0 or n0 == 0:
        return float("nan")
    rank = 1
    rank_sum = 0.0
    i = 0
    while i < len(scores):
        j = i
        while j < len(scores) and scores[j][0] == scores[i][0]:
            j += 1
        avg_rank = (rank + rank + (j - i) - 1) / 2.0
        rank_sum += avg_rank * sum(label for _, label in scores[i:j])
        rank += j - i
        i = j
    return (rank_sum - n1 * (n1 + 1) / 2.0) / (n1 * n0)


def summarize(values):
    if not values:
        return {}
    vals = sorted(values)
    return {
        "min": vals[0],
        "p10": vals[min(len(vals) - 1, int(len(vals) * 0.10))],
        "median": statistics.median(vals),
        "p90": vals[min(len(vals) - 1, int(len(vals) * 0.90))],
        "max": vals[-1],
        "mean": statistics.mean(vals),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pairs_tsv", nargs="?", default="isabl-test.pairs.tsv")
    ap.add_argument("--column", default="concordance")
    ap.add_argument("--limit", type=int, default=12)
    args = ap.parse_args()

    rows = []
    with Path(args.pairs_tsv).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        sample_a_col = "#sample_a" if "#sample_a" in reader.fieldnames else "sample_a"
        for row in reader:
            try:
                expected = float(row["expected_relatedness"])
                score = float(row[args.column])
            except Exception:
                continue
            if expected not in (-1.0, 1.0):
                continue
            rows.append((expected, score, row[sample_a_col], row["sample_b"]))

    pos = [score for expected, score, _, _ in rows if expected == 1.0]
    neg = [score for expected, score, _, _ in rows if expected == -1.0]
    roc = auc([(score, 1 if expected == 1.0 else 0) for expected, score, _, _ in rows])

    print(f"column={args.column}")
    print(f"rows={len(rows)} pos={len(pos)} neg={len(neg)}")
    print(f"auroc={roc:.6f}")
    print(f"pos<0.99={sum(v < 0.99 for v in pos)}")
    print(f"pos<0.9={sum(v < 0.9 for v in pos)}")
    print(f"pos<0.8={sum(v < 0.8 for v in pos)}")
    print(f"pos<0.7={sum(v < 0.7 for v in pos)}")
    print(f"pos<0.6={sum(v < 0.6 for v in pos)}")
    print(f"neg<0.4={sum(v < 0.4 for v in neg)}")
    print(f"neg>0.4={sum(v > 0.4 for v in neg)}")
    print(f"pos_summary={summarize(pos)}")
    print(f"neg_summary={summarize(neg)}")
    print("worst_pos=")
    for item in sorted((score, a, b) for expected, score, a, b in rows if expected == 1.0)[:args.limit]:
        print(item)
    print("highest_neg=")
    for item in sorted(((score, a, b) for expected, score, a, b in rows if expected == -1.0), reverse=True)[:args.limit]:
        print(item)


if __name__ == "__main__":
    main()
