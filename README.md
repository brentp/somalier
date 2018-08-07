# somalier: checking cancer sample matching from BAMs/CRAMs

Existing software for checking relatedness requires jointly-called germ-line variant calls,
but cancer projects have only somatic calls between tumor-normal pairs.

`somalier` makes checking any number of samples for identity easy:

```
somalier -o results/my-project --sites $vcf *.cram
```

where sites is a VCF of variant sites (provided by somalier for hg19).
The interactive output that's produced (similar to [peddy](https://github.com/brentp/peddy))
makes it fast and easy to detect mismatched samples and sample-swaps.

An optional argument lets the user specify expected groups of samples.


## Usage

```
somalier [options] <bam/cram>...

Arguments:
  <bam/cram> file(s) for samples of interest.

Options:

  -s --sites <vcf>        vcf file with lines of sites to use for relatedness estimation.
  -t --threads <int>      optional number of processors to use for parallelization.
  -f --fasta <reference>  path to reference fasta file.
  -g --groups <path>      optional path to expected groups of samples (e.g. tumor normal pairs).
                          specified as comma-separated groups per line e.g.:
                            normal1,tumor1a,tumor1b
                            normal2,tumor2a
  -o --output <prefix>    output prefix for results.
```

## How it works

`somalier` takes a list of known sites. Even a few hundred (or dozen) sites can be a very
good indicator of relatedness. The best sites are those with a population allele frequency
close to 0.5 as that maximizes the probability that any 2 samples will differ.

In order to quickly calculate genotypes at these sites, `somalier` assays the exact base
[without using pileup](https://brentp.github.io/post/no-pile/). It also parallelizes across
sites with as many threads as requested. In addition, it uses [hts-nim](https://github.com/brentp/hts-nim)
which is a very fast wrapper of [htslib](https://htslib.org).

For each sample-pair, it reports:
1. IBS0 -- the number of sites where one sample is hom-ref and another is hom-alt
2. IBS2 -- the number of sites where the samples have the same genotype
3. shared-hets -- the number of sites where both samples are heterozygotes
4. shared-hom-alts -- the number of sites where both samples are homozygous alternate

These are used to calculate [relatedness](https://en.wikipedia.org/wiki/Coefficient_of_relationship)
and a measure of relatedness that is unaffected by loss-of-heterozygosity that is often seen in some 
cancers. The interactive output allows toggling between any of these measures.


## Example

![example](readme-fig.png)


## Other Work

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5499645/

https://academic.oup.com/bioinformatics/article/33/4/596/2624551


## TODO

1. link plots
2. guess which samples are normal.
[x] label as related/unrelated instead of 0 or 1
