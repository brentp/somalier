f=supplemental-file-3.txt
echo -e "n\ttp\tfp\tfn\tstrict" > $f
./subset_sites --ped 3tissue-v2.ped          -n 10 -n 20 -n 40 -n 100 -n 200 -n 400 -n 1000 -n 2000 -n 4000 -n 8000 -n 16000 --vcf subset-sites/sites.hg38.rna.v2.vcf.gz subset-sites/somalier-v2/*.somalier >> $f
./subset_sites --strict --ped 3tissue-v2.ped -n 10 -n 20 -n 40 -n 100 -n 200 -n 400 -n 1000 -n 2000 -n 4000 -n 8000 -n 16000 --vcf subset-sites/sites.hg38.rna.v2.vcf.gz subset-sites/somalier-v2/*.somalier >> $f

python plot-subset-sites.py $f
