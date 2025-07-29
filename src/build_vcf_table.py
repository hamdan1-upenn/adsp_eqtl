import argparse
import sys
import pandas as pd
from helper_functions import *
import os

parser = argparse.ArgumentParser(description='Build a table from VCF file used for eQTL analysis.')
parser.add_argument('-v','--vcf-file', help='Path to the VCF file.')
parser.add_argument('-a','--rnaseq-file', help='Path to the RNA-seq anndata file. Anndata object must have ENSEMBL gene IDs as the var index.')
parser.add_argument('-o', '--outfile', help='Output file for tabularized VCF File (don\'t include .gz suffix. Program will gzip file post-processing).')
parser.add_argument('-r', '--rosmap', action='store_true', help='Use this if the VCF file is from ROSMAP. This will automatically remove the unnecessary suffix from sample IDs')
working_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
bin_path = os.path.join(working_dir, 'bin')

args = parser.parse_args()

cmd=f'bash {bin_path}/tabularize_vcf.sh {bin_path}/bcftools-1.21/bcftools {args.vcf_file} {args.outfile}'
os.system(cmd)


if args.rosmap:
    snp_df = pd.read_csv(f'{args.outfile}.gz', sep='\t')
    snp_df.columns = [col.split('-B')[0] for col in snp_df.columns]
    snp_df.columns = [col.split('-U')[0] for col in snp_df.columns]
    snp_df.to_csv(f'{args.outfile}.gz', sep='\t', index=False, compression='gzip')
