dev
===
+ if a sample had > 1 allele that was neither REF nor ALT at a given site, it was assigned
  an `unknown` genotype. This was too stringent for deep sequencing so it was changed to a
  proportion (> 0.04 [or 1 in 25 alleles]) #7

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
