# Cancer, Concordance, Contamination

Version 0.3.2 and higher of somalier estimate contamination in several ways.

## Background

Especially in cancer genomics, it's important to estimate inter-sample contamination of a normal into a tumor sample since the normal is used to call somatic variants. There are tools like `conpair` and `verifyBamId` that do this well. 

In cancer samples with structural deletions and duplications, the allele balance (AB) is not expected to be 0.5 for non-homozygous sites. Conpair and CHARR therefore look only at sites likely to be homozygous (AB close to 0 or 1). A simple method that also uses this fact has also been present in `somalier` as p_middling_ab for many years.

It also now supports `concordance` in a way that is robust to copy-number variation. That will be output by the `relate` sub-command for all pairs of samples. Tumor-normal pairs should have a concordance > 0.6 (should be > 0.95 for high-quality samples) and unrelated or non-self pairs should have a concordance that approaches 0, but below 0.4.

## Contamination

This was present in somalier as the proportion of sites with an allele balance (AB) in the ranges (0, 0.1) or (0.9, 1.0), which are considered nearly homozygous. This works well because it avoids looking at heterozygotes which are skewed by copy-number variation and loss-of-heterozygosity--an observation utilized by CHARR and CONPAIR. This is still available as `p_middling_ab` per-sample.

As of 0.3.2, there is also a CHARR(-like) estimate of contamination. This works well and uses population allele frequencies (which are available in the somalier sites files) to scale the probability that reads with the other allele in nearly homozygous sites are from contamination. This is a good single-sample estimate until contamination levels get too high. It appears in the html output.

Somalier also introduces a `contamination` sub-command that implements a CONPAIR-like estimate for a tumor-normal pair. Like conpair, and like CHARR, this focuses on apparent homozygous sites. Our conpair-like implementation then looks at those sites in the tumor. Since the sites are common germ-line variants, it is expected that their appearance at any level indicates contamination from a different person in the tumor sample.

The contamination sub-command can be run all-vs-all, but neither it nor conpair are designed for that use. It is recommended to run it with a tumor-normal pair using the `--tumor-normal-pair` (`-p`) flag. It will do all vs all without that flag if given `*.somalier`, but that will require scrutinizing along with other `somalier` output to diagnose.


**NOTE**: for cancer samples, set `-d 20` or higher to help with the contamination estimates. The default is lower for germ-line samples.

For cancer cohorts, it's useful to look at the `relate` html output and:
- set the Y-axis of the relatedness plot to `Concordance`
- click the `Contamination QC` preset at the bottom to set the axis of the sample plot on the right.
