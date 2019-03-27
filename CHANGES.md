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
