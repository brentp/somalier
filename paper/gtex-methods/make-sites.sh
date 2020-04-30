#wget https://storage.googleapis.com/gtex_analysis_v8/reference/gencode.v26.GRCh38.genes.gtf
#wget https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_exon_reads.parquet
PATH=$PATH:~/Projects/src/bedtools2/bin
mkdir -p tmp
export TMPDIR=$(pwd)/tmp/
set -euo pipefail

D=../simons-dn/
ped=$D/10065_66.ped
vcf=$D/simons_p231.hg38.bcf
ref=~/Data/GRCh38_full_analysis_set_plus_decoy_hla.fa

python median-by-exon.py
python to-bed.py gencode.v26.GRCh38.genes.gtf  gtex.exon-counts.quantiles.tsv > gtex.v8.exon-counts.quantiles.bed

# $12 is 40th percentile
awk '$12 > 10' gtex.v8.exon-counts.quantiles.bed | cut -f 1-3 \
	| sort -k1,1 -k2,2n \
	| bedtools merge > gtex.v8.40pct.gt10.bed

pslivar expr \
	--fasta $ref \
	--ped $ped \
	-v $vcf \
	-x /uufs/chpc.utah.edu/common/HIPAA/u6000771/Data/LCR-hs38.bed.gz \
	--pass-only \
	--info "variant.ALT.length == 1 && INFO.AC > 6 && variant.REF.length == 1 && variant.ALT[0].length == 1 && variant.ALT[0] != '*'" \
	--trio "snp_dn:(kid.GQ > 10 && dad.GQ > 10 && mom.GQ > 10 && kid.DP > 8 && dad.DP > 8 && mom.DP > 8 && kid.AB > 0.1 && kid.AB < 0.9) && ((kid.het && mom.hom_ref && dad.hom_ref) || (!kid.het && mom.hom_ref && dad.hom_alt) || (!kid.het && mom.hom_alt && dad.hom_ref))" \
	  | bcftools view --threads 8 -O b -o simons-ac.gt2.bcf

bcftools index --threads 12 simons-ac.gt2.bcf                                                                                                         
slivar make-gnotate --prefix dn --field AC:dn_AC simons-ac.gt2.bcf 

somalier_dev find-sites \
    -i ~/Projects/2020/gtex-somalier/gtex.v8.40pct.gt10.bed \
    -x ~/Data/LCR-hs38.bed.gz \
    --min-AN 20000 ~/Data/gnomad/3/gnomad.hg38.no-vep.AFgt0.01.vcf.gz \
    --snp-dist 6000 --gnotate-exclude \
    dn.zip \
    && mv sites.vcf.gz sites.hg38.rna.vcf.gz

