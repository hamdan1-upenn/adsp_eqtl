import argparse
import sys
from helper_functions import *
import pdb

parser = argparse.ArgumentParser(add_help=False,description='Run eQTL analysis using memento. use -h for particular formatting detail')
parser.add_argument('-s','--snp-file', help='Path to the SNP file.This is the tabularized VCF file made with build_vcf_table.py.')
parser.add_argument('-a','--rnaseq-file', help='Path to the RNA-seq anndata file. Anndata object must have ENSEMBL gene IDs as the var index.')
parser.add_argument('-v','--covariate-file',help='path to CSV file containing covariate dataframe. ALL VALUES IN THESE COLUMNS MUST BE NUMERIC.')
parser.add_argument('-k', '--subject-key', help='Column name in AnnData obs dataframe that corresponds to Subject IDs in the SNP file')
parser.add_argument('-o', '--output-file', help='Output file for eQTL results.')
parser.add_argument('--celltype-column', help='Column name in AnnData obs dataframe that corresponds to cell type')
parser.add_argument('-c', '--celltype', help='Cell type to be subsetted for eQTL analysis.')
parser.add_argument('--snp-ids', help='EITHER a comma-separated list of SNP ID with no spaces OR a path to a file with SNP_IDs on each line to be used for eQTL analysis. `Use chr_pos_ref_alt` format. CANNOT use this and --snp-range in the same call')
parser.add_argument('--snp-range', help='A range within the genome where all SNPs between the range will be collected for analysis. Use chr_start_end format. CANNOT use this and --snp-ids in the same call')
parser.add_argument('--gene-range', help='A range within the genome where all genes between the range will be collected for analysis. Use chr_start_end format')
parser.add_argument('-h','--help', action='store_true', help='Show this help message and exit')

args= parser.parse_args()

if args.snp_ids and args.snp_range:
    print('ERROR: cannot use --snp-ids and --snp-range in the same call. Exiting...')
    sys.exit(1)

if args.help:
    print_help_msg(parser)

check_file_exists(args.snp_file)
check_file_exists(args.rnaseq_file)
check_file_exists(args.covariate_file)

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
gtf_path = os.path.join(working_dir, 'bin', 'Homo_sapiens.GRCh38.113.gtf')
bin_path = os.path.join(working_dir, 'bin')


snp_df = read_snp_into_df(args.snp_file)
try:
    snp_df['snp_id'] = snp_df['CHROM'].astype(str) + '_' + snp_df['POS'].astype(str) + '_' + snp_df['REF'].astype(str) + '_' + snp_df['ALT'].astype(str)
except:
    print('ERROR: could not build  snp_id column for snp dataframe. Check the format of the tabularized vcf file (use -h for more info)\nIf VCF file was originally used, tabularized file can be found in the vcf file\'s directory within \'tmp\'')

# load in snp IDs
if args.snp_ids:
    if os.path.isfile(args.snp_ids):
        with open(args.snp_ids, 'r') as file:
            target_snps = [line.strip() for line in file.readlines()]
    else:
        target_snps = args.snp_ids.split(',')
elif args.snp_range:
    chrom,start,end = args.snp_range.split('_')
    snp_df = snp_df[(snp_df.CHROM==chrom) & (snp_df.POS >= int(start)) & (snp_df.POS <= int(end))]
    if snp_df.empty:
        print('ERROR: no SNPs found in the specified range. Exiting...')
        sys.exit(1)
    target_snps = snp_df.snp_id.tolist()
else:
    print('ERROR: no target SNPs specified. Exiting...')
    sys.exit(1)

# find snps found in snp dataframe and exit if no target snps match
if not set(target_snps) & set(snp_df.snp_id):
    print('ERROR: no matching SNP IDs found in the SNP dataframe. Exiting...')
    sys.exit(1)
else:
    snp_list = list(set(target_snps) & set(snp_df.snp_id))
    print(f'Found {len(snp_list)} matching SNP IDs in the SNP dataframe')

# read in anndata object
try:
    print('READING IN RNA-SEQ DATA...')
    adata = sc.read_h5ad(args.rnaseq_file)
    print('FINISHED READING RNA-SEQ DATA')
except:
    print('ERROR: could not read the rnaseq file. Check the format of the h5ad file')

# ensure celltype and subject key columns are in the adata object
check_column_name(adata.obs, args.subject_key)
check_column_name(adata.obs, args.celltype_column)

subject_key = args.subject_key
celltype_column = args.celltype_column

#create covariate dataframe
cov_df = pd.read_csv(args.covariate_file)
# check if subject key is in the covariate dataframe
if not (set(cov_df.iloc[:,0]) & set(snp_df.columns)):
    print('ERROR: first column in covariate dataframe does not match the subject IDs in the vcf dataframe. Exiting...')
    sys.exit(1)
cov_df.index = cov_df.iloc[:,0]
cov_df = cov_df.iloc[:,1:]

snp_df.index = snp_df.snp_id
snp_df = snp_df[list(set(snp_df.columns) & set(cov_df.index) & set(adata.obs[subject_key].unique()))]
if snp_df.empty:
    print('ERROR: no matching Sample IDs found in the SNP dataframe. Exiting...')
    sys.exit(1)

# check if celltype is in the adata celltype column
if args.celltype not in adata.obs[celltype_column].unique():
    print('ERROR: cell type not found in the adata object.')
    print(f'cell types in column \"{celltype_column}\":\n{adata.obs[celltype_column].unique()}\nExiting...')
    sys.exit(1)
else:
    print(f'Found {args.celltype} in the adata object')
    adata = adata[adata.obs[celltype_column] == args.celltype]

# check if subject IDs for VCF file are in anndata object
if not (set(adata.obs.subject_key.unique()) & set(snp_df.columns)):
    print('ERROR: VCF file and subject key column have no matching sample IDs. Check the format of the VCF file and the subject key IDs')
    sys.exit(1)

print('BUILDING GTF REFERENCE DATAFRAME...')
gtf = build_gtf_ref(gtf_path)
print('FINISHED BUILDING GTF REFERENCE DATAFRAME')


gene_snp_pair = []
if args.snp_range:
    if args.gene_range:
        chrom, start, end = args.gene_range.split('_')
    genes_in_range = find_snp_genes_within_range(gtf,chrom, start, end)
    for snp in snp_list:
        if len(genes_in_range) > 0:
            genes = [[snp, gene] for gene in genes_in_range]
            gene_snp_pair.extend(genes)
else:# create snp pair dataframe for each gene within 500kb of the snp
    for snp in snp_list:
        genes_within_500kb = list(find_snp_genes_within500kb(gtf, snp))
        if len(genes_within_500kb) > 0:
            genes = [[snp, gene] for gene in genes_within_500kb]
            gene_snp_pair.extend(genes)

subject_keys = adata.obs[subject_key].unique()

gene_snp_pair = [item for item in gene_snp_pair if item is not None]
gene_snp_df = pd.DataFrame(gene_snp_pair, columns=['Position', 'Gene'])
gene_snp_df = gene_snp_df.drop_duplicates().reset_index(drop=True)

gene_snp_df = gene_snp_df.rename(columns={'Position':'SNP','Gene':'gene'})
gene_snp_df = gene_snp_df[['gene','SNP']]


dfs = []
for snp in target_snps:
    tmp_snp_df = snp_df.iloc[snp_df.index == snp]
    if tmp_snp_df.empty:
        print(f'SNP {snp} not found in the SNP dataframe. continueing...')
        continue
    # remove columns with NaN values
    tmp_snp_df = tmp_snp_df.dropna(axis=1)
    overlap_ids = set(adata.obs[subject_key].drop_duplicates().tolist()) & set(tmp_snp_df.columns) & set(cov_df.index)  
    adata = adata[adata.obs[subject_key].isin(list(overlap_ids))]
    #adata = adata[:, adata.var.index.isin(gene_snp_df['gene'].unique())]

    memento_ready_snp_df = tmp_snp_df.loc[:,list(overlap_ids)]
    cov = cov_df.loc[list(overlap_ids),:]
    
    memento_ready_gene_snp_df = gene_snp_df[gene_snp_df['gene'].isin(list(adata.var.index))]
    memento_ready_gene_snp_df = memento_ready_gene_snp_df[memento_ready_gene_snp_df['SNP'].isin(list(memento_ready_snp_df.index))]
    if (adata.obs.empty | memento_ready_gene_snp_df.empty | memento_ready_snp_df.empty):
        print('ERROR: one or more of the dataframes are empty. Exiting...')
        print(f'adata: {adata.obs.empty}, memento_ready_gene_snp_df: {memento_ready_gene_snp_df.empty}, memento_ready_snp_df: {memento_ready_snp_df.empty}')
        sys.exit(1)
    try:
        eqtl_results = memento.run_eqtl(adata=adata,
                        snps=memento_ready_snp_df.T,
                        cov=cov,
                        gene_snp_pairs=memento_ready_gene_snp_df,
                        donor_column=subject_key,
                        num_cpu=32,
                        num_blocks=1
                        )
        dfs.append(eqtl_results)

    except Exception as e:
        print(f'Error running eQTL analysis for SNP {snp}: {e}')
        continue

eqtl_results = pd.concat(dfs, axis=0)
eqtl_results.index = eqtl_results.gene.astype(str)
eqtl_results = pd.merge(eqtl_results, gtf[['gene_name']], left_index=True, right_index=True, how='left')
eqtl_results.to_csv(args.output_file)


#eqtl_results = memento.run_eqtl(adata=adata, snps=memento_ready_snp_df.T, cov=cov, gene_snp_pairs=memento_ready_gene_snp_df, donor_column=subject_key, num_cpu=1, num_blocks=1 )

#snp_fpath = '/project/guomiclab/hamdanz_projects/adsp_r4_rosmap/data/adsp_r4_rosmap/chr19.rosmap.r4.wgs.biallelic.genotypes.mac10.corrected.txt.gz'

#snp_df = pd.read_csv(snp_fpath, sep='\t')
#snp_df['snp_id'] = snp_df['CHROM'].astype(str) + '_' + snp_df['POS'].astype(str) + '_' + snp_df['REF'].astype(str) + '_' + snp_df['ALT'].astype(str)

#adata_fpath = '/project/guomiclab/hamdanz_projects/adsp_r4_rosmap/data/rnaseq/rosmap_raw_rnaseq_labeled.h5ad'
#adata.var.index = adata.var.gene_name
#adata = adata[:, ~adata.var.index.duplicated(keep='last')]





