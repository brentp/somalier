#./src/somalier find-sites /data/gemini_install/data/gemini_data/gnomad.exomes.r2.1.tidy.bcf
set -euo pipefail

if [[ -f b37tohg19.chain ]]; then
	echo OK
else
wget -q https://github.com/broadgsa/gatk/raw/master/public/chainFiles/b37tohg19.chain
fi

picard LiftoverVcf \
I=sites.vcf.gz \
O=sites.hg19.tmp.vcf.gz \
CHAIN=b37tohg19.chain \
REJECT=rejected_variants.37to19.vcf \
R=/data/human/hg19.fa

if [[  -f hg19ToHg38.over.chain.gz ]]; then
	echo "OK"
else
wget -q http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
fi

picard LiftoverVcf \
I=sites.hg19.tmp.vcf.gz \
O=sites.hg38.vcf.gz \
RECOVER_SWAPPED_REF_ALT=true \
CHAIN=hg19ToHg38.over.chain.gz \
REJECT=rejected_variants.19to38.vcf \
R=/data/human/hg38.fa


if [[ -f hg38ToHg19.over.chain.gz ]]; then
	echo "OK"
else
wget -q http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
fi

picard LiftoverVcf \
I=sites.hg38.vcf.gz \
O=sites.hg19.vcf.gz \
RECOVER_SWAPPED_REF_ALT=true \
CHAIN=hg38ToHg19.over.chain.gz \
REJECT=rejected_variants.38to19.vcf \
R=/data/human/hg19.fa


zcat sites.hg19.vcf.gz | grep -Pv "samples of|non_cancer|controls|non_neuro" | perl -pe 's/^chr//' | bgzip -c > sites.GRCh37.vcf.gz
zcat sites.hg38.vcf.gz | grep -Pv "samples of|non_cancer|controls|non_neuro" | perl -pe 's/^chr//' | bgzip -c > sites.hg38.nochr.vcf.gz

zgrep -cv ^# sites.GRCh37.vcf.gz sites.hg38.nochr.vcf.gz sites.hg38.vcf.gz
