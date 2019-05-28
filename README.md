# somalier: extract informative sites and evaluate relatedness and quality-control

Existing software for checking relatedness requires jointly-called germ-line variant calls,
but cancer projects have only somatic calls between tumor-normal pairs.

`somalier` makes checking any number of samples for identity easy **directly from the alignments**:

The first step is to extract sites. This is parallelizable by sample (here using `xargs`, but this can also be done via cluster or cloud):
```
ls *.cram | xargs -I{} -P 12 "somalier extract -d extracted/ --sites sites.vcf.gz -f /data/human/g1k_v37_decoy.fa {}"
```

`--sites` is a VCF of known polymorphic sites in VCF format. A good set is provided in
the [releases](https://github.com/brentp/somalier/releases) but any set of common variants will work.


The next step is to calculate relatedness on the extracted data:

```
somalier relate --ped $pedigree extracted/*.somalier
```
This will create text and interactive HTML output that's produced (similar to [peddy](https://github.com/brentp/peddy))
makes it fast and easy to detect mismatched samples and sample-swaps.

Note that the `somalier relate` command runs extremely quickly (10 seconds for 600 samples) so it's possible
to add/remove samples or adjust a pedigree file and re-run iteratively.

For example to add the **n + 1th** samples, just run `somalier extract` on the new sample and then re-use
the already extracted data from the `n` original samples.


## Install

get a static binary from [here](https://github.com/brentp/somalier/releases)

Users can also get a docker image [here](https://hub.docker.com/r/brentp/somalier/tags)
which contains htslib and a somalier binary ready-for-use.

## How it works

`somalier` takes a list of known polymorphic sites. Even a few hundred (or dozen) sites
can be a very good indicator of relatedness. The best sites are those with a population
allele frequency close to 0.5 as that maximizes the probability that any 2 samples will differ.
A list of such sites is provided in the [releases](https://github.com/brentp/somalier/releases)
for GRCh37 and hg38.

In order to quickly calculate genotypes at these sites, `somalier` assays the exact base.
The extraction step is done directly from the bam/cram files 1 sample at a time.

The `rel` step is run on the output of the `extract` commands. It runs extremely quickly
so that new samples can be added and compared.

For each sample-pair, it reports:
1. IBS0 -- the number of sites where one sample is hom-ref and another is hom-alt
2. IBS2 -- the number of sites where the samples have the same genotype
3. shared-hets -- the number of sites where both samples are heterozygotes
4. shared-hom-alts -- the number of sites where both samples are homozygous alternate

These are used to calculate [relatedness](https://en.wikipedia.org/wiki/Coefficient_of_relationship)
and a measure of relatedness that is unaffected by loss-of-heterozygosity that is often seen in some 
cancers. The interactive output allows toggling between any of these measures.

It also reports depth information and the count of `HET`, `HOM_REF`, `HOM_ALT`, and `unknown` genotypes for each sample
along with a number of metrics that are useful for general QC.

## Example

![example](https://user-images.githubusercontent.com/1739/43783575-4863f13c-9a1f-11e8-9cf8-622f784edc69.png)

Here, each point is a pair of samples. We can see that the expected identical sample-pairs (e.g. tumor-normal pairs) specified by the user
and drawn in red mostly cluster together on the right. Unrelateds cluster on the lower left. The sample-swaps are the blue points that cluster with
the red. In the somalier output, the user can **hover to see which sample-pairs are involved each point**

## License

`somalier` is free and unrestricted for non-commercial use. For commercial use, please contact [bpedersen@base2genomics.com]

## Other Work

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5499645/

https://academic.oup.com/bioinformatics/article/33/4/596/2624551


## Acknowledgement

This work was motivated by interaction and discussions with Preeti Aahir and several
early users who provided valuable feedback.
