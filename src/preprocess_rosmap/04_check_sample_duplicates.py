import scanpy as ac
import pandas as pd
from collections import Counter
import os


snp_fpath = '/project/guomiclab/hamdanz_projects/adsp_r4_rosmap/data/adsp_r4_rosmap/chr19.rosmap.r4.wgs.biallelic.genotypes.mac10.txt.gz'

snp_df = pd.read_csv(snp_fpath, sep='\t')

rna_id_match = [col.split('-B')[0] for col in snp_df.columns]
rna_id_match = [col.split('-U')[0] for col in rna_id_match]

dupes = [item for item, count in Counter(rna_id_match).items() if count > 1]
full_dupe_id = [dupe for dupe in snp_df.columns if dupes[0] in dupe]











