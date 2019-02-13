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
