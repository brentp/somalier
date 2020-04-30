import pandas as pd
import numpy as np

df = pd.read_parquet("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_exon_reads.parquet")
#, engine='fastparquet')
print(df.shape)
print(df.head())

df.quantile(np.arange(0.05, 0.951, 0.05), axis=1).T.to_csv('gtex.exon-counts.quantiles.tsv',
        sep="\t", float_format="%.2f", index_label="exon")

