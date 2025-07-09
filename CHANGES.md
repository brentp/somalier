v0.3.0
======
+ fix compilation with nim > 2 and some sex inference stuff

v0.2.19
=======
+ [relate/infer] fix check that would prevent some inference (#123)

v0.2.18
=======
+ [find_sites] handle empty alternate alleles (#121)
+ [find_sites] add --output-vcf option 

v0.2.17
=======
+ allow setting bam sample name via ENV when RG is missing (#115)

v0.2.16
=======
+ make find-sites more efficient and more robust (#96)

v0.2.14
=======
+ minor fixes to find-sites
+ allow setting env variable `SOMALIER_REPORT_ALL_PAIRS` to force reporting of all sample-pairs (#76)
+ improve readme (via @zztin)

v0.2.13
=======
+ add "Heterozygosity rate" as a per-sample metric to the html output. (Thanks Irenaeus and Kelly for the suggestion)
+ fix inference for some cases. obvious parent-child pairs were sometimes missed.

v0.2.12
=======
+ add checkbox to HTML to scale IBS0, IBS2, etc by number of sites shared by the samples. this almost always
  results in a scaling that is better across (pairs of) samples.
+ ancestry: allow globs for ancestry files (#59)

v0.2.11
=======
+ more informative error message on bad sample name (#53)
+ allow setting SOMALIER_AB_HOM_CUTOFF to change which calls are considered hom-ref (#56)
+ adjust (fix) relatedness calculation which was off when the number of shared sites was low (#55).
  many thanks to @fgvieira who found and diagnosed this problem. This change adds a `hets_ab` column
  which is the count of times sample `a` was het and `b` was **not** unknown + the times sample `b` was het
  and `a` was not unknown; this is mostly not needed except to (re) calculate relatedness but is reported
  in the text output for completeness.

v0.2.10
=======
+ fix extra output column in pairs.tsv (#47)
+ change output file name of ancestry to include "somalier"
+ fix for gvcf with empty alts (#46)
+ add include regions and exclude sites to find-sites
+ add --min-ab option to somalier relate to limit het sites to `min_ab..(1-min_ab)`. default is 0.3
+ html output for sample plot defaults to number of het sites on X (was hom-alt)
+ better estimates in `somalier ancestry` when incoming samples are different ancestry from training (thousand genomes samples)

v0.2.9
======
+ support multi-sample GVCF and fix some GVCF cases (thanks @ameynert for implementing)
+ also fixes some edge-cases with GVCFs

v0.2.8
======
+ html output has a list of pre-sets to auto-select informative X, Y axes for the sample plot
+ add --infer flag to somalier relate to allow inferring relatedness.
  this accompanies a change in the .samples.tsv output so that it can be used as a pedigree file
+ add --sample-prefix option to extract and corresponding (multi-)option to relate. So, given
  a cohort with DNA and RNA where samples have identical IDs (SM tags) in the DNA and RNA, can use
  somalier as:
  ```
  somalier extract -d DNA --sample-prefix DNA- ...
  somalier extract -d RNA --sample-prefix RNA- ...

  somalier relate --sample-prefix DNA- --sample-prefix RNA- DNA/*.somalier RNA/*.somalier ...
  ```

  and it will show the samples that have
  matching IDs after stripping the prefixes as "identical".


v0.2.7
======
+ new subcommand `ancestry` to predict ancestry using a simple neural network on the somalier
  sketches. creates an interactive html output and a text file
+ fix for "Argument list too long" on huge cohorts (#37)
+ sub-sample .pairs.tsv output for huge cohorts -- only for unrelated samples.
+ better sub-sampling of html output

v0.2.6
======
+ for large cohorts (>1K samples) the html output is now usable.
  it randomly subsets samples that should be and are unrelated.
+ better error messages for bad input
+ inspect environment variable: `SOMALIER_ALLOWED_FILTERS` so that users can
  give a comma-delimited list of FILTERs that should be allowed (by default only PASS and RefCall
  variants are considered. This is useful for some GVCF formats. 
  
v0.2.5
======
+ handle more types of GVCF (#27, thanks @holtjma)
+ handle VCFs without depth (AD) information. this enables extracting VCFs with only genotypes such as
  files converted from array information (#31, thanks @asazonov)

v0.2.4
======
+ unify genotyping between all code-paths (thanks Filipe)
+ if both groups and pedigree information are specified, they correctly share information (#26)
+ relax allele balance to hom-ref is < 0.04 and hom-alt > 0.96 (was 0.02 and 0.98 respectively).
+ support for GVCF (#27)

v0.2.3
======
+ calculate relatedness correctly for samples with parent-ids specified
  when the parents are not actually in the pedigree file.
+ use bit-vectors to calculate relatedness. this gives about a 250X speedup.
  with this code, I can now evaluate relatedness for 3756 in under 30 seconds on my laptop.
+ better scaling of X and Y depth
+ use final RG as the sample id in relate
+ output expected relatedness in .pairs.tsv file
+ fix ref/alt (a/b-allele ordering for VCF) this was a bug that caused problems when comparing
  samples extracted from VCF files to other samples extracted from BAM/CRAM files. Thanks very 
  much to Filipe and Sergio for finding this issue and providing several test-cases. (if you
  have previously downloaded the thousand genomes files from zenodo, please update to the latest).

v0.2.2
======
+ add a default output directory
+ static build with libcurl

v0.2.1 
======
+ fix hover in html
+ add --unknown flag for `somalier relate` to set unknown genotypes to hom-ref (useful when merging single-sample VCFs).
+ change sites to be alphabetical by allele so that they are the same between genome builds
+ add version to .somalier files created with extract -- these will not be compatible with those made with v0.2.0. I don't
  forsee a backwards incompatible change like this one in the near future.
+ sites files for hg38 and GRCh37 are compatible. That is, we can extract sites from bams or vcfs from samples aligned to GRCh37
  reference and accurately calculate relatedness on files extracted from samples aligned to hg38.
+ better HTML performance for large numbers of samples by sub-sampling individiuals that are expected to be unrelated and that 
  have a calculated relatedness < 0.09.
+ add a `depthview` sub-command to plot the depth of each sample along each chromosome.
+ much nicer html and several fixes thanks to Joe Brown

v0.2.0
======
This was a large re-write of `somalier`. The command-line usage is backwards incompatible (but
should not change moving forward). There is now a per-sample extract step:
```
somalier extract -d extracted/ -s $sites_vcf -f $fasta $sample.cram
```

followed by a relate step:

```
somalier relate --ped $ped extracted/*.somalier
```

This enables parallelization by sample across nodes and the resulting, extracted, binary "somalier"
files are only ~240KB per sample so reading them is nearly instant and the `relate` step
runs in 10 seconds for my 603-sample test-case which makes adjusting pedigree files or removing samples
and re-running a much faster process.

somalier extract can also take a (multi-sample) VCF and create an idential "somalier" file
for cases when a VCF is available. 

v0.1.5
======
+ add experimental contamination estimate. this simply prints to stderr the sample and
  inferred source (another sample) of contamination along with the estimated level of 
  contamination and the number of sites used to estimate it.
+ fix threading bug with large numbers of samples. 
+ more lenient ped file parsing ("Female" will be recognized in sex column and
  "Affected" in phenotype column).
+ the html output now allows selecting a single sample to be highlighted in the plot
  this allows finding a sample of interest in a large cohort.
+ the output now includes a new metric for proportion of sites with an allele balance
  > 0.02 and < 0.2 or > 0.8 and < 0.98. this turns out to be a nice QC (high is bad)
+ for low coverage or targetted sites, sometimes `nan` values would stop the entire 
  html page from working; this has been fixed.
+ make sure all reported relationships are plotted in correct colors (#14)
+ plotting fixes (#15)

v0.1.4
======
+ if a file ending with ".list" is given as an argument (instead of .bam, .cram), it can contain
  paths to the alignment files and optionally the indexes. e.g.
  ```
https://abc/path/to/aaa.bam https://abc/indexes/path/aaa.bam.bai
https://abc/path/to/bbb.bam https://abc/indexes/path/bbb.bam.bai
```
  These can be space, comma, or tab-delimited.

v0.1.3
======
+ if a sample had > 1 allele that was neither REF nor ALT at a given site, it was assigned
  an `unknown` genotype. This was too stringent for deep sequencing so it was changed to a
  proportion (> 0.04 [or 1 in 25 alleles]) #7
+ for samples with sparse coverage, e.g. from targetted sequencing projects, mean depth is
  not very informative because it gets washed out by all the zero-depth sites. The new columns:
  `gt_depth_mean`, `gt_depth_std`, gt_depth_skew` report the values for the depth at genotyped
  sites--those meeting the depth requirement (default of 7).

v0.1.2
======
+ allow lower-case reference alleles in case of masked genomes (see #5)
+ set relatedness values < -1.5 to -1.5 in the plot
+ fix bug that affected relatedness calcs especially in RNA-Seq
+ add more diagnostic values (allele-balance and number of non ref/alt bases)

v0.1.1
======
+ fix bug in plot labels
+ better inter-plot interaction in html
