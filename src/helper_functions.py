import os
import sys
import pandas as pd

#check if file exists. print status
def check_file_exists(file_path):
    if os.path.isfile(file_path):
        print(f'{file_path} exists! continuing...')
        return
    else:
        print(f'{file_path} does not exists! exiting...')
        sys.exit(1)


#prints help message from parser
def print_help_msg(parser):
    format_instructions = """
    Tabularized SNP File Column Format:

    CHROM\tPOS\tID\tREF\tALT\tRSID\tRSID\tAC\tAF\tSAMPLE_ID1\tSAMPLE_ID2\t...\tSAMPLE_IDn
    
    -- sample ID values should be 0, 1 or 2.
    -- If inputting a VCF file directly, use the --vcf flag. This will convert the VCF file into a tabularized format.
    
    """
    print(format_instructions)
    parser.print_help()
    sys.exit(0)


#read tabularized snp file into a dataframe
def read_snp_into_df(snp_fpath):
    try:
        snp_df = pd.read_csv(snp_fpath,sep='\t')
        return snp_df
    except:
        print('Error reading tabularized SNP file to dataframe')


