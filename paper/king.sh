time plink2 --make-bed --vcf <(bcftools view --threads 4 1kg-ccdg.bcf) --out plink-1kg
time king --kinship -b plink-1kg.bed

time somalier extract -s ~/Data/sites.hg38.vcf.gz -d 1kg-somalier/ -f $ref 1kg-ccdg.bcf
time somalier relate  -o 1kg 1kg-somalier/*.somalier

python paper/king-fig.py paper/1kg.pairs.tsv paper/king.kin0
