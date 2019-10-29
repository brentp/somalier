# somalier: extract informative sites, evaluate relatedness, and perform quality-control on BAM/CRAM/BCF/VCF/GVCF

[![Actions Status](https://github.com/brentp/somalier/workflows/Docker%20Image%20CI/badge.svg)](https://github.com/brentp/somalier/actions)

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
This will create text and interactive HTML output (similar to [peddy](https://github.com/brentp/peddy))
that makes it fast and easy to detect mismatched samples and sample-swaps.

Note that the `somalier relate` command runs extremely quickly (< 2 seconds for 600 samples and ~1 minute for 4,500 samples) so it's possible
to add/remove samples or adjust a pedigree file and re-run iteratively.

For example to add the **n + 1th** samples, just run `somalier extract` on the new sample and then re-use
the already extracted data from the `n` original samples.

## Usage

The usage is also described above. Briefly, run:
```
somalier extract -d cohort/ --sites sites.hg38.vcf.gz -f $reference $sample.bam
```
for each sample to create a small binary file of the ref and alt counts for the variants listed
in sites.hg38.vcf.gz.

for a vcf, run:
```
somalier extract -d cohort/ --sites sites.hg38.vcf.gz -f $reference $cohort.bcf
```

Then run:
```
somalier relate --ped $pedigree_file cohort/*.somalier
```
This will create an html file for QC in a few seconds. 

Note that if a new sample is added to the cohort, it's only necessary to perform
the `extract` step on that sample and then run the (fast) `relate` step again with all
of the extracted files.

## VCF

`somalier` can `extract` from a multi or single-sample VCF or a GVCF. This will be much faster, in cases where it's available,
this would look like:

```
somalier extract -d extracted/ --sites sites.vcf.gz -f /data/human/g1k_v37_decoy.fa joint.vcf.gz
```

following this, there will be a `$sample.somalier` file for each sample in the `joint.vcf.gz`

Note that `somalier` uses the `AD` field to extract depth information. If that FORMAT field is not present in the
header, then it will use the genotypes only and use a total depth of 20 (10,10 for heterozygote), etc.


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

The `relate` step is run on the output of the `extract` commands. It runs extremely quickly
so that new samples can be added and compared. It uses 3 bit-vectors per sample for hom-ref,
het, hom-alt. Each bitvector is a sequence of 64 bit integers where each bit is set if
the variant at that index in the sample is for example, heterozygous. With this setup,
we can use fast bitwise operations and [popcount](https://en.wikichip.org/wiki/population_count)
hardware instructions to calculate relatedness extremely quickly.

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


## Ancestry Estimate

You can use a background set of extracted `*.somalier` files, for example from thousand genomes along with labels for
those samples to predict ancestry for a given set (or single sample). This would look like:

```
python scripts/ancestry-predict.py --labels scripts/ancestry-labels-1kg.tsv --samples $MY_SAMPLES/*.somalier --backgrounds 1kg-somalier/*.somalier > sample-ancestries.txt
```

The ancestry-predict.py and the `scripts/ancestry-labels-1kg.tsv` are in the somalier repository.
The somalier files for thousand genomes can be downloaded from [here](https://zenodo.org/record/3479773/files/1kg.somalier.tar.gz?download=1)
These were created from the thousand genomes high coverage data from [here](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20190425_NYGC_GATK/)
Note that these will work for either GRCh37 or hg38 as long as you use the most recent sites files distributed with `somalier`.


## Usage

Usage is intentionally very simple and running `somalier extract` or `somalier relate` will give sufficient help for nearly
all cases.

By default `somalier` will only consider variants that have a "PASS" or "RefCall" FILTER. To extend this list, set
the environment variable `SOMALIER_ALLOWED_FILTERS` to a comma-delimited list of additional filters to allow.
 
## License

`somalier` is free and unrestricted for non-commercial use. For commercial use, please contact [bpedersen@base2genomics.com]

## Other Work

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5499645/

https://academic.oup.com/bioinformatics/article/33/4/596/2624551


## Acknowledgement

This work was motivated by interaction and discussions with Preeti Aahir and several
early users who provided valuable feedback.
