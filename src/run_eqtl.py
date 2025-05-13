import argparse
import sys
from helper_functions import *

parser = argparse.ArgumentParser(add_help=False,description='Run eQTL analysis using memento. use -h for particular formatting detail')
parser.add_argument('-s','--snp-file', help='Path to the SNP file. This can be in either VCF or tabular format (see -h for formatting). If VCF, use the --vcf flag')
parser.add_argument('--vcf', action='store_true', help='Use this if SNP file is in VCF format. This will be use to convert VCF into tabularized DataFrame.')
parser.add_argument('-a','--rnaseq-file', help='Path to the RNA-seq anndata file')
parser.add_argument('-k', '--subject-key', help='Column name in AnnData obs dataframe that corresponds to Subject IDs in the SNP file')
parser.add_argument('-o', '--output-dir', help='Path to the output directory')
parser.add_argument('--celltype-column', help='Column name in AnnData obs dataframe that corresponds to cell type')
parser.add_argument('-c', '--celltype', help='Cell type to be subsetted for eQTL analysis. If not specified, results will be generated for all cell types')
parser.add_argument('--snp-ids', help='EITHER acomma-separated list of SNP ID with no spaces OR a path to a file with SNP_IDs on each line to be used for eQTL analysis. `Use chr_pos_ref_alt` format')
parser.add_argument('-h','--help', action='store_true', help='Show this help message and exit')


args= parser.parse_args()

if args.help:
    print_help_msg(parser)

check_file_exists(args.snp_file)
check_file_exists(args.rnaseq_file)

print('Importing Packages...')
import scanpy as sc
import numpy as np
import pandas as pd
import memento
import pdb 
import os
import shutil
import subprocess
print('Packages finished importing')

working_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
gtf_path = os.path.join(working_dir, 'data','metadata', 'Homo_sapiens.GRCh38.113.gtf')
bin_path = os.path.join(working_dir, 'bin')

print('READING VCF/SNP DATA...')
if args.vcf:
    print('CONVERTING VCF TO TABULARIZED FORMAT...')

    #remove any existing tmp directory
    if os.path.exists(f'{os.path.dirname(args.snp_file)}/tmp'):
        shutil.rmtree(f'{os.path.dirname(args.snp_file)}/tmp')

    tabularize_args=['bash', f'{bin_path}/tabularize_vcf.sh', f'{bin_path}/bcftools', args.snp_file]

    result=subprocess.run(tabularize_args, capture_output=True, text=True)
    
    if result.returncode != 0:
        print('Error Tabularizing VCF. File must be in standard VCF format')
        sys.exit(result.returncode)
    else:
        snp_fpath = f'{os.path.dirname(args.snp_file)}/tmp/tabularized_vcf.tsv'
        snp_df = read_snp_into_df(snp_fpath)
        shutil.rmtree(f'{os.path.dirname(args.snp_file)}/tmp')

else:
    snp_df = read_snp_into_df(args.snp_file, sep='\t')

try:
    snp_df['snp_id'] = snp_df['CHROM'].astype(str) + '_' + snp_df['POS'].astype(str) + '_' + snp_df['REF'].astype(str) + '_' + snp_df['ALT'].astype(str)
except:
    print('ERROR: could not build  snp_id column for snp dataframe. Check the format of the tabularized vcf file (use -h for more info)\nIf VCF file was originally used, tabularized file can be found in the vcf file\'s directory within \'tmp\'')



#snp_fpath = '/project/guomiclab/hamdanz_projects/adsp_r4_rosmap/data/adsp_r4_rosmap/chr19.rosmap.r4.wgs.biallelic.genotypes.mac10.corrected.txt.gz'

if os.path.isfile(args.snp_ids):
    with open(args.snp_ids, 'r') as file:
        target_snps = [line.strip() for line in file.readLines()]
else:
    target_snps = args.snp_ids.split(',')


#snp_df = pd.read_csv(snp_fpath, sep='\t')
#snp_df['snp_id'] = snp_df['CHROM'].astype(str) + '_' + snp_df['POS'].astype(str) + '_' + snp_df['REF'].astype(str) + '_' + snp_df['ALT'].astype(str)

adata_fpath = '/project/guomiclab/hamdanz_projects/adsp_r4_rosmap/data/rnaseq/rosmap_raw_rnaseq_labeled.h5ad'
adata = sc.read_h5ad(adata_fpath)

sc_batch_snp_batch_column_ids = ('batch','subject_key')

GTF_PATH = '/project/guomiclab/hamdanz_projects/adsp_r4_rosmap/data/metadata/Homo_sapiens.GRCh38.113.gtf'
gtf = pd.read_csv(GTF_PATH, sep='\t', header= None, comment = '#')

gtf = gtf[gtf[2]=='gene']
gtf = gtf.rename(columns={3:'start',4:'end'})
gtf['ensembl_id'] = gtf[8].str.extract('gene_id "([^"]+)"').values.astype(str)
gtf['gene_name'] = gtf[8].str.extract('gene_name "([^"]+)"').values.astype(str)
gtf['chrom'] = gtf[0].values.astype(str)
gtf = gtf[['chrom','start', 'end', 'ensembl_id','gene_name']]
gtf=gtf[gtf.notnull().all(1)]
gtf = gtf[gtf['gene_name'] !='nan']
gtf.index = list(gtf['ensembl_id'])


adata.var.index = adata.var.gene_name
adata = adata[:, ~adata.var.index.duplicated(keep='last')]

#create a gene_snp dataframe pairing for every gene within 500kb of an snp
def find_snp_genes_within500kb(df, snp_id):
    df['chrom']=df['chrom'].replace('chr','')
    chrom, pos = snp_id.split('_')[0:2]
    chrom=chrom.replace('chr','')
    pos=int(pos)
    start = max(pos-500000,1)
    end = min(pos+500000,max(df['start']))
    df = df[df['chrom'] == chrom]
    gene_ref = df[(df['start'] >= start) & (df['start'] <= end)]
    return(gene_ref['gene_name'].values)

snp_list = target_snp
gene_snp_pair = []
for snp in snp_list:
    genes = list(find_snp_genes_within500kb(gtf, snp))
    if len(genes) > 0:
        genes = [[snp, gene] for gene in genes]
        gene_snp_pair.extend(genes)

subject_keys = adata.obs['subject_key'].unique()
snp_df.index = snp_df.snp_id

gene_snp_pair = [item for item in gene_snp_pair if item is not None]
gene_snp_df = pd.DataFrame(gene_snp_pair, columns=['Position', 'Gene'])
gene_snp_df = gene_snp_df.drop_duplicates().reset_index(drop=True)

#create covariate dataframe for memento
cov_df = adata.obs[['subject_key', 'msex', 'age_death']]
cov_df = cov_df.drop_duplicates().reset_index(drop=True)
cov_df.index = cov_df.subject_key
cov_df.drop(columns='subject_key',inplace=True)

#convert all 90+ to just 90 for covariate df
def cap_age_at_90(age):
    if age == '90+':
        age='90'
    return age

cov_df['age_death'] = cov_df['age_death'].map(lambda x: cap_age_at_90(x))

gene_snp_df = gene_snp_df.rename(columns={'Position':'SNP','Gene':'gene'})
gene_snp_df = gene_snp_df[['gene','SNP']]



for ct in adata.obs.major_cell_type.unique():
    ct_adata = adata[adata.obs.major_cell_type==ct]
    overlap_ids = set(ct_adata.obs.subject_key.drop_duplicates().tolist()) & set(snp_df.columns) & set(snp_df.columns) & set(cov_df.index)  
    ct_adata = ct_adata[ct_adata.obs['subject_key'].isin(list(overlap_ids))]
    memento_ready_snp_df = snp_df.loc[:,list(overlap_ids)]
    memento_ready_snp_df = memento_ready_snp_df.fillna(0)
    cov = cov_df.loc[list(overlap_ids),:]
    memento_ready_gene_snp_df = gene_snp_df[gene_snp_df['gene'].isin(list(ct_adata.var.index))]
    eqtl_results = memento.run_eqtl(adata=ct_adata,
                        snps=memento_ready_snp_df.T,
                        cov=cov,
                        gene_snp_pairs=memento_ready_gene_snp_df,
                        donor_column='subject_key',
                        num_cpu=1,
                        num_blocks=1
                        )
    eqtl_results.to_csv(f'/project/guomiclab/hamdanz_projects/adsp_r4_rosmap/output/{ct}_chr19_rosmap_r4_adsp.csv')


