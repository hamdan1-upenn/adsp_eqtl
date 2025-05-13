import scanpy as ac
import pandas as pd
import os

snp_fpath = '/project/guomiclab/hamdanz_projects/adsp_r4_rosmap/data/adsp_r4_rosmap/chr19.rosmap.r4.wgs.biallelic.genotypes.mac10.txt.gz'
snp_outpath = '/project/guomiclab/hamdanz_projects/adsp_r4_rosmap/data/adsp_r4_rosmap/chr19.rosmap.r4.wgs.biallelic.genotypes.mac10.corrected.txt.gz'

snp_df = pd.read_csv(snp_fpath, sep='\t')

#there are apendices to the sample names that do not match the adata sample match key that need to be removed
snp_df.columns = [col.split('-B')[0] for col in snp_df.columns]
snp_df.columns = [col.split('-U')[0] for col in snp_df.columns]

snp_df.to_csv(snp_outpath, sep='\t', index=False, compression='gzip')

adata_fpath = '/project/guomiclab/hamdanz_projects/adsp_r4_rosmap/data/rnaseq/rosmap_raw_rnaseq_labeled.h5ad'
adata = sc.read_h5ad(adata_fpath)

def cap_age_at_90(age):
    if age == '90+':
        age='90'
    return age

cov['age_death'] = cov['age_death'].map(lambda x: cap_age_at_90(x))

adata.var.index = adata.var.gene_name
adata = adata[:, ~adata.var.index.duplicated(keep='last')]

adata.write_h5ad('/project/guomiclab/hamdanz_projects/adsp_r4_rosmap/data/rnaseq/rosmap_raw_rnaseq_labeled_memento_corrected.h5ad')




