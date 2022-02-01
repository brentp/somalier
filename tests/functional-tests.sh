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

prefix_fn() {
    set -ex
    $exe extract --sample-prefix AA -s tests/test_sites.vcf -f tests/test.fa -d test_prefix_A tests/gt_only.vcf.gz
    $exe extract --sample-prefix BB -s tests/test_sites.vcf -f tests/test.fa -d test_prefix_B tests/gt_only.vcf.gz
    SOMALIER_SAMPLE_RATE=0.0 $exe relate -o test_prefix_A/out --sample-prefix BB --sample-prefix AA test_prefix_B/*.somalier test_prefix_A/*.somalier
}

export -f prefix_fn

run check_prefix prefix_fn
assert_exit_code 0

# checks that final column, expected relatedness is 1.
assert_equal "1" $(awk 'NR == 1 && NR == 1.0' test_prefix_A/out.pairs.tsv | wc -l) 
rm -rf test_prefix_A test_prefix_B

assert_in_stderr "got sampling rate 0.0 from env var"

run check_gvcf_no_alt $exe extract -s tests/x.gvcf.gz -f tests/test.fa tests/x.gvcf.gz
assert_exit_code 0

