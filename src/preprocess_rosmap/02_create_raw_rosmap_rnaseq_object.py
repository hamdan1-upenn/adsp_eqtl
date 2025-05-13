import scanpy as sc
import pandas as pd
import anndata as ad
from concurrent.futures import ProcessPoolExecutor
import os
import shutil

INPUT_DIR = '/project/guomiclab/hamdanz_projects/LOY_integrated_sea-ad_rosmap/data/rosmap/snrnaseq/adata'

OUTPUT_FILE= '/project/guomiclab/hamdanz_projects/adsp_r4_rosmap/data/rnaseq/raw_rosmap_rnaseq.h5ad'

GTF_PATH = '/project/guomiclab/hamdanz_projects/adsp_r4_rosmap/data/metadata/Homo_sapiens.GRCh38.113.gtf'

PT_META_PATH = '/project/guomiclab/hamdanz_projects/adsp_r4_rosmap/data/metadata/ROSMAP_ADSPID_merged_metadata.csv'
RNASEQ_IND_META = '/project/guomiclab/hamdanz_projects/adsp_r4_rosmap/data/metadata/ROSMAP_rnaseq_metadata.csv'

h5ad_files = [os.path.join(INPUT_DIR, f) for f in os.listdir(INPUT_DIR) if f.endswith(".h5ad")]

temp_dir='/project/guomiclab/hamdanz_projects/adsp_r4_rosmap/data/tmp_h5ad'

# in order to concat on disk, must be in csr format:
def convert_to_csr(file):
    adata = ad.read_h5ad(file)
    if not adata.X.format == "csr":  # Check if the matrix is not in CSR format
        print(f"Converting {file} to CSR format...")
        adata.X = adata.X.tocsr()  # Convert to CSR format
        temp_file = os.path.join(temp_dir, os.path.basename(file))
        adata.write_h5ad(temp_file)  # Write to a temporary file
        return temp_file
    else:
        print(f"{file} is already in CSR format.")
        return file 

with ProcessPoolExecutor() as executor:
    results = list(executor.map(convert_to_csr, h5ad_files))

tmp_files = [os.path.join(temp_dir, os.path.basename(file)) for file in h5ad_files]
ad.experimental.concat_on_disk(tmp_files, OUTPUT_FILE, join='outer')

adata = ad.read_h5ad(OUTPUT_FILE)

gtf = pd.read_csv(GTF_PATH, sep='\t', header= None, comment = '#')

gtf = gtf[gtf[2]=='gene']
gtf['ensembl_id'] = gtf[8].str.extract('gene_id "([^"]+)"').values.astype(str)
gtf['gene_name'] = gtf[8].str.extract('gene_name "([^"]+)"').values.astype(str)
gtf['chrom'] = gtf[0].values.astype(str)
gtf = gtf[['chrom','ensembl_id','gene_name']]
gtf=gtf[gtf.notnull().all(1)]
gtf = gtf[gtf['gene_name'] !='nan']
gtf.index = list(gtf['ensembl_id'])

adata.var['ensembl_id'] = list(adata.var.index)
adata.var= pd.merge(adata.var, gtf[['gene_name','chrom']], left_index=True, right_index=True, how='left')

rnaseq_id_meta = pd.read_csv(RNASEQ_IND_META)
rnaseq_id_meta['cellID'] = rnaseq_id_meta.libraryBatch.astype(str) + '#' + rnaseq_id_meta.cellBarcode.astype(str)


pt_meta = pd.read_csv(PT_META_PATH)

obs_meta = pd.merge(rnaseq_id_meta, pt_meta, left_on='individualID', right_on='individualID', how='left')
obs_meta.index = obs_meta['cellID']

adata.obs = pd.merge(adata.obs, obs_meta, left_index=True, right_index=True, how='left')

adata = adata[~(adata.obs.subject_key.isnull())]

adata.write_h5ad(OUTPUT_FILE)


