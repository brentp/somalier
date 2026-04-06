# Concordance Experiments

## Baseline

- Baseline metric: `inferred_hom_concordance` alone.
- Evaluation command:
  - `bash eval.sh`
  - `python3 scripts/evaluate_concordance.py isabl-test.pairs.tsv`
- Baseline summary on `isabl-test.pairs.tsv` before this round of changes:
  - `AUROC`: `0.999889`
  - `expected_relatedness == 1.0`: mean `0.9956`, median `1.0`
  - `expected_relatedness == -1.0`: mean `0.5427`, median `0.5395`
  - `neg < 0.4`: `0 / 3622`
- Conclusion: ranking was already excellent, but absolute calibration was poor for unrelated pairs.

## Offline Experiments

- Pure multiplicative hom-alt penalty:
  - Formula family: `concordance * hom_alt_raw_clipped^a`
  - Effect: pushed almost all unrelated pairs below `0.4`, but dropped many true `1.0` pairs too far from `1`.
- Pure thresholded hom-alt penalty:
  - Formula family: `concordance - w * max(0, t - hom_alt_raw_clipped)`
  - Best simple version kept positive mean high (`~0.987`) and pushed `3615 / 3622` negatives below `0.4`.
  - Failure mode: several true DNA/RNA matches with low hom-alt signal were penalized too hard.
- `p_middling_ab` as a global scale:
  - Formula family: multiply the score by `(1 - p_middling_ab)` or use it as a direct scale.
  - Effect: did not help enough by itself; it mostly shrank all scores.
- `p_middling_ab` as a penalty amplifier:
  - Formula family: `concordance - w * max(0, t - hom_alt_raw_clipped) * (1 + b * mean_p_middling_ab)`
  - Better than a global scale. It kept the positive median close to `1` while making low hom-alt/noisy pairs fall faster.
- Small relatedness regularizer:
  - Added `- d * (1 - max(relatedness, 0))` to catch the last clean false positives where both concordance-style signals stayed high.
  - This improved the unrelated tail without changing the top of the true-match distribution much.
- Information-scaled penalties:
  - Formula family: multiply the penalty by a function of `n`, for example `sqrt(min(1, n / N0))` or `sqrt(n / (n + N0))`.
  - Effect: this recovered several low-information DNA/RNA true matches, but it also reopened too many unrelated RNA-heavy pairs above `0.4`.
  - Conclusion: not worth the tradeoff on this dataset.
- Upper-range stretch above the discordance cutoff:
  - Formula family: leave scores `<= 0.4` unchanged, and remap the interval `(0.4, 1.0]` with a concave transform.
  - Best simple version used a square-root stretch anchored at `0.4`.
  - Effect: improved calibration for true matches without changing how many unrelated pairs stayed below `0.4`.

## Six More Experiments

- Experiment 1: `base` only.
  - Formula: `concordance = inferred_hom_concordance`
  - Result: `pos < 0.9 = 1`, `neg > 0.4 = 3622`, `pos < 0.8 = 0`
  - Conclusion: far too permissive for unrelated pairs.
- Experiment 2: `base` plus upper-range stretch only.
  - Formula: stretch `base` above `0.4`
  - Result: `pos < 0.9 = 0`, `neg > 0.4 = 3622`, `pos < 0.8 = 0`
  - Conclusion: even worse for negatives; stretch alone cannot calibrate discordance.
- Experiment 3: `min(base, hom_alt)` plus stretch.
  - Formula: `stretch(min(base, hom_alt))`
  - Result: `pos < 0.9 = 27`, `neg > 0.4 = 5`, `pos < 0.8 = 6`
  - Conclusion: simple, but over-penalizes true matches.
- Experiment 4: thresholded hom-alt penalty without `p_middling_ab` or relatedness.
  - Formula: `stretch(base - max(0, 0.72 - hom_alt))`
  - Result: `pos < 0.9 = 4`, `neg > 0.4 = 37`, `pos < 0.8 = 0`
  - Conclusion: the shape for positives is good, but it needs a noise-aware penalty amplifier.
- Experiment 5: thresholded hom-alt penalty with `p_middling_ab`, no relatedness.
  - Formula family: `stretch(base - max(0, t - hom_alt) * (1 + b * pm))`
  - Best result in this family: `t = 0.72`, `b = 12`, `stretch exponent = 0.4`
  - Result: `pos < 0.9 = 4`, `neg > 0.4 = 4`, `pos < 0.8 = 3`
  - Conclusion: this matches the current target with one fewer moving part.
- Experiment 6: thresholded hom-alt penalty with relatedness but no `p_middling_ab`.
  - Formula family: `stretch(base - max(0, t - hom_alt) - d * (1 - relatedness))`
  - Best result in this family: `neg > 0.4 = 10`, `pos < 0.9 = 4`, `pos < 0.8 = 3`
  - Conclusion: relatedness alone is less useful than `p_middling_ab` for separating the hard negatives.

## Five More Experiments

- Experiment 7: simple threshold `0.70`, amplifier `10`, square-root stretch.
  - Formula: `stretch(base - max(0, 0.70 - hom_alt) * (1 + 10 * pm))`, with `stretch(x) = anchor + span * sqrt(...)`
  - Result: `pos < 0.9 = 5`, `neg > 0.4 = 6`, `pos < 0.8 = 3`, `pos < 0.7 = 2`
  - Conclusion: too weak on both fronts.
- Experiment 8: simple threshold `0.75`, amplifier `10`, square-root stretch.
  - Result: `pos < 0.9 = 6`, `neg > 0.4 = 4`, `pos < 0.8 = 4`, `pos < 0.7 = 3`
  - Conclusion: negatives are fine, but positives degrade too much.
- Experiment 9: simple threshold `0.70`, amplifier `10`, cube-root stretch.
  - Result: `pos < 0.9 = 4`, `neg > 0.4 = 6`, `pos < 0.8 = 2`, `pos < 0.7 = 1`
  - Conclusion: good positive shape, but not enough negative suppression.
- Experiment 10: simple threshold `0.75`, amplifier `10`, cube-root stretch.
  - Result: `pos < 0.9 = 4`, `neg > 0.4 = 4`, `pos < 0.8 = 3`, `pos < 0.7 = 2`
  - Conclusion: this is the cleanest simple-number formula that still matches the target.
- Experiment 11: simple threshold `0.70`, amplifier `15`, square-root stretch.
  - Result: `pos < 0.9 = 5`, `neg > 0.4 = 4`, `pos < 0.8 = 3`, `pos < 0.7 = 3`
  - Conclusion: stronger amplifier fixes negatives, but harms positives more than Experiment 10.

## Five More Experiments Again

- Experiment 12: `base` only, no penalty and no stretch.
  - Result: `pos < 0.6 = 0`, `neg > 0.4 = 3622`
  - Conclusion: unusable because every unrelated pair stays high.
- Experiment 13: threshold `0.75`, amplifier `10`, no stretch.
  - Result: `pos < 0.6 = 3`, `neg > 0.4 = 4`, `pos < 0.7 = 4`, `pos < 0.8 = 6`
  - Conclusion: the simplest viable metric so far, but it loses one extra hard positive relative to the current version.
- Experiment 14: threshold `0.70`, amplifier `10`, no stretch.
  - Result: `pos < 0.6 = 2`, `neg > 0.4 = 6`, `pos < 0.7 = 3`, `pos < 0.8 = 4`
  - Conclusion: too permissive for negatives.
- Experiment 15: threshold `0.75`, amplifier `10`, square-root stretch.
  - Result: `pos < 0.6 = 2`, `neg > 0.4 = 4`, `pos < 0.7 = 3`, `pos < 0.8 = 4`
  - Conclusion: a good simple runner-up, but still worse than cube root on the hardest positives.
- Experiment 16: `min(base, hom_alt)`, no stretch.
  - Result: `pos < 0.6 = 5`, `neg > 0.4 = 5`, `pos < 0.7 = 12`, `pos < 0.8 = 32`
  - Conclusion: too aggressive and not competitive.

## Ten More Experiments

- Experiment 17: threshold `0.70`, amplifier `5`, no stretch.
  - Result: `pos < 0.6 = 2`, `neg > 0.4 = 7`, `pos < 0.7 = 3`
  - Conclusion: not enough negative suppression.
- Experiment 18: threshold `0.75`, amplifier `10`, no stretch.
  - Result: `pos < 0.6 = 3`, `neg > 0.4 = 4`, `pos < 0.7 = 4`
  - Conclusion: viable, but loses extra positives.
- Experiment 19: threshold `0.70`, amplifier `15`, square-root stretch.
  - Result: `pos < 0.6 = 2`, `neg > 0.4 = 4`, `pos < 0.7 = 3`
  - Conclusion: a decent simple fallback, but still weaker than cube root.
- Experiment 20: threshold `0.75`, amplifier `10`, square-root stretch.
  - Result: `pos < 0.6 = 2`, `neg > 0.4 = 4`, `pos < 0.7 = 3`
  - Conclusion: also viable, but not the best positive recovery.
- Experiment 21: threshold `0.70`, amplifier `5`, cube-root stretch.
  - Result: `pos < 0.6 = 0`, `neg > 0.4 = 7`, `pos < 0.7 = 0`
  - Conclusion: best positive recovery, but too many negatives remain above `0.4`.
- Experiment 22: threshold `0.70`, amplifier `10`, cube-root stretch.
  - Result: `pos < 0.6 = 1`, `neg > 0.4 = 6`, `pos < 0.7 = 1`
  - Conclusion: close, but still worse than the best balanced candidates.
- Experiment 23: threshold `0.70`, amplifier `15`, cube-root stretch.
  - Grid result: `pos < 0.6 = 1`, `neg > 0.4 = 4`, `pos < 0.7 = 2`
  - Real eval result: `pos < 0.6 = 1`, `neg > 0.4 = 4`, `pos < 0.7 = 2`
  - Conclusion: best simple-number candidate so far.
- Experiment 24: threshold `0.75`, amplifier `5`, cube-root stretch.
  - Result: `pos < 0.6 = 0`, `neg > 0.4 = 7`, `pos < 0.7 = 1`
  - Conclusion: same issue as Experiment 21; too permissive for negatives.
- Experiment 25: threshold `0.70`, amplifier `15`, simple linear stretch with slope `1.5`.
  - Result: `pos < 0.6 = 2`, `neg > 0.4 = 4`, `pos < 0.7 = 3`
  - Conclusion: understandable, but weaker than cube root.
- Experiment 26: threshold `0.75`, amplifier `10`, simple linear stretch with slope `1.5`.
  - Result: `pos < 0.6 = 3`, `neg > 0.4 = 4`, `pos < 0.7 = 3`
  - Conclusion: not competitive with the better simple candidates.

## Simpler Penalty

- Additive penalty instead of multiplicative penalty:
  - Formula: `penalty = max(0, 0.70 - hom_alt + 2 * p_middling_ab)`
  - Real eval result: `pos < 0.6 = 0`, `pos < 0.7 = 0`, `pos < 0.8 = 2`, `neg > 0.4 = 6`
  - Conclusion: this is meaningfully simpler and still acceptable under the relaxed budget, even though it allows two more negatives above `0.4`.

## Implemented Candidate

- Current implemented formula:
  - `base = inferred_hom_concordance`
  - `hom_alt = clip(raw_hom_alt_concordance, 0, 1)`
  - `pm = mean(sample_a.p_middling_ab, sample_b.p_middling_ab)`
  - `penalized = clip(base - max(0, 0.70 - hom_alt + 2 * pm), 0, 1)`
  - `concordance = penalized` for `penalized <= 0.4`
  - `concordance = 0.4 + 0.6 * cbrt((penalized - 0.4) / 0.6)` for `penalized > 0.4`
- Reason:
  - Preserves the current interpretation of the low-concordance region.
  - Uses `p_middling_ab` as a knob only when hom-alt concordance is suspiciously low.
  - Removes the relatedness term entirely, which makes the metric simpler.
  - Uses simple, defensible constants: `0.70`, `2`, and cube root.
  - Lifts the compressed upper range so true matches land closer to `1.0`.

## Current Result

- Evaluation command:
  - `bash eval.sh`
  - `python3 scripts/evaluate_concordance.py isabl-test.pairs.tsv`
- Result after the anchored upper-range stretch:
  - `AUROC`: `0.999819`
  - `expected_relatedness == 1.0`: mean `0.9906`, median `1.0`
  - `expected_relatedness == -1.0`: mean `0.00328`, median `0.0`
  - `pos < 0.99`: `15 / 119`
  - `pos < 0.9`: `4 / 119`
  - `pos < 0.8`: `2 / 119`
  - `pos < 0.7`: `0 / 119`
  - `pos < 0.6`: `0 / 119`
  - `neg < 0.4`: `3616 / 3622`
  - `neg > 0.4`: `6 / 3622`
- Comparison to the previous implemented version:
  - `neg > 0.4` increased from `4` to `6`
  - `pos < 0.7` improved from `2` to `0`
  - `pos < 0.6` improved from `1` to `0`
  - This is a deliberate tradeoff in favor of a much simpler penalty.

## Notes

- `eval.sh` needed `DYLD_FALLBACK_LIBRARY_PATH=/opt/homebrew/lib` in this environment so the binary could load `libhts.dylib`.
- The evaluation target here is calibration, not just ranking. AUROC stayed high for many formulas, so the useful metric was mainly:
  - how close the `1.0` group stayed to `1`
  - how many `-1.0` rows dropped below `0.4`
