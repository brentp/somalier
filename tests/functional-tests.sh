#!/bin/bash
test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest

. ssshtest

set -o nounset

nim c -d:debug  -d:useSysAssert -d:useGcAssert --lineDir:on --debuginfo --boundChecks:on -x:on src/somalier
set +e
exe=./src/somalier

run check_help_works $exe --help
assert_exit_code 0
assert_in_stdout "Commands:"

#wget -O - https://github.com/brentp/somalier/files/3412455/sites.GRCh37.vcf.gz | zcat - | head -200 > tests/test_sites.vcf
# make a GT only VCF

set -e
awk 'BEGIN{FS=OFS="\t"} ($0 ~ /^#CHROM/) { print $0"\tFORMAT\ttest_sample"; next } NR == 3 { print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"} $0 ~ /^#/ { print $0 } $0 !~ /^#/ { if(NR!=180) {print $0"\tGT\t0/1"} else { print $0"\tGT\t./." }}' \
	tests/test_sites.vcf | bgzip -c > tests/gt_only.vcf.gz
tabix tests/gt_only.vcf.gz
set +e


run check_gt_only_extract $exe extract -s tests/test_sites.vcf -f tests/test.fa -d test_gt_only tests/gt_only.vcf.gz
assert_in_stderr "[somalier] FORMAT field 'AD' not found for depth information. using genotype only"
assert_exit_code 0
assert_in_stderr "[somalier] found 49 sites"

rm -rf test_gt_only
#rm -f tests/gt_only.vcf.gz
