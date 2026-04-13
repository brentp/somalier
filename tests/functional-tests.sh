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

prefix_fn() {
	pre
    $exe extract --sample-prefix AA -s tests/test_sites.vcf -f tests/test.fa -d test_prefix_A tests/gt_only.vcf.gz
    $exe extract --sample-prefix BB -s tests/test_sites.vcf -f tests/test.fa -d test_prefix_B tests/gt_only.vcf.gz
    $exe relate -s tests/test_sites.vcf -o test_prefix_A/out --sample-prefix BB --sample-prefix AA test_prefix_B/*.somalier test_prefix_A/*.somalier
}

export -f prefix_fn

run check_prefix prefix_fn
assert_exit_code 0

# checks that final column, expected relatedness is 1.
assert_equal "1" $(awk 'NR == 1 && NR == 1.0' test_prefix_A/out.pairs.tsv | wc -l) 
assert_equal "1" $(head -n 1 test_prefix_A/out.samples.tsv | grep -c "contamination_charr")
rm -rf test_prefix_A test_prefix_B


run check_gvcf_no_alt $exe extract -s tests/x.gvcf.gz -f tests/test.fa tests/x.gvcf.gz
assert_exit_code 0

prefix_contam_fn() {
    mkdir -p test_contam
    $exe extract --sample-prefix A- -s tests/test_sites.vcf -f tests/test.fa -d test_contam_A tests/gt_only.vcf.gz
    $exe extract --sample-prefix B- -s tests/test_sites.vcf -f tests/test.fa -d test_contam_B tests/gt_only.vcf.gz
    $exe contamination -s tests/test_sites.vcf test_contam_A/*.somalier test_contam_B/*.somalier -o test_contam/out.tsv
}

bad_sites_contam_fn() {
    mkdir -p test_contam
    $exe extract --sample-prefix BAD- -s tests/test_sites.vcf -f tests/test.fa -d test_contam_bad tests/gt_only.vcf.gz
    awk 'BEGIN{n=0} /^#/ { print; next } n < 5 { print; n++ }' tests/test_sites.vcf > tests/bad_contam_sites.vcf
    $exe contamination -s tests/bad_contam_sites.vcf test_contam_bad/*.somalier -o test_contam/bad.tsv
}

export -f prefix_contam_fn
export -f bad_sites_contam_fn

run check_contamination_sample_rows prefix_contam_fn
assert_exit_code 0
assert_in_stderr "wrote test_contam/out.samples.tsv and test_contam/out.pairs.tsv"
assert_equal "3" $(wc -l < test_contam/out.samples.tsv | tr -d ' ')
assert_equal "3" $(wc -l < test_contam/out.pairs.tsv | tr -d ' ')
assert_equal "#sample_name\tanchor_sample\tn_sites_usable\tcontamination_mle" "$(head -n 1 test_contam/out.pairs.tsv)"
assert_equal "1" $(awk '$1 == "A-test_sample" && $2 == "B-test_sample"' test_contam/out.pairs.tsv | wc -l | tr -d ' ')
assert_equal "1" $(awk '$1 == "B-test_sample" && $2 == "A-test_sample"' test_contam/out.pairs.tsv | wc -l | tr -d ' ')
rm -rf test_contam_A test_contam_B test_contam

run check_contamination_bad_sites bad_sites_contam_fn
assert_exit_code 1
assert_in_stderr "expected 5 from --sites"
rm -rf test_contam_bad test_contam
rm -f tests/bad_contam_sites.vcf
