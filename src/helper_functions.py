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

#check if file exists. print status
def check_dir_exists(dir_path):
    if os.path.isdir(dir_path):
        print(f'{dir_path} exists! continuing...')
        return
    else:
        print(f'{dir_path} does not exists! exiting...')
        sys.exit(1)


#prints help message from parser
def print_help_msg(parser):
    format_instructions = """
    Tabularized SNP File Column Format:

    CHROM\tPOS\tID\tREF\tALT\tRSID\tRSID\tAC\tAF\tSAMPLE_ID1\tSAMPLE_ID2\t...\tSAMPLE_IDn
    
    -- sample ID values should be 0, 1 or 2.
    
    INPUTING TARGET SNPs:

    -- using `--snp-range` will require the input chr_start_end format. for example, chr13_133000_233000 will get every snp between chr13:133000 and chr13:233000.
    -- using `--snp-ids` will require the input snp_id format. for example, chr13_133000_A_T will get the snp with id chr13_133000_A_T.
        -- when using `--snp-ids`, you can input a comma separated list with no spaces (ex. chr13_133000_A_T,chr13_133001_A_T,...) or a file path that has one snp_id per line.
    
    
    COVARIATE FILE FORMAT:

    subject_key\tcov1\tcov2\t...\tcovn

    -- subject_key should match the IDs in the SNP dataframe
    -- all covariates should be numeric
    -- subject key must be the first column
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

##check if column name exists in dataframe
def check_column_name(df, col_name):
    if col_name not in df.columns:
        print(f'Column name {col_name} not found in dataframe. Exiting...')
        sys.exit(1)
    else:
        print(f'Column name {col_name} found in dataframe. Continuing...')


#build gtf reference dataframe
def build_gtf_ref(GTF_PATH):
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
    return gtf


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
    return(gene_ref['ensembl_id'].values)

#create a gene_snp dataframe pairing for every gene within 500kb of an snp
def find_snp_genes_within_range(df,chrom, start, end):
    chrom=chrom.replace('chr','')
    start=int(start)
    end=int(end)
    df = df[df['chrom'] == chrom]
    gene_ref = df[(df['start'] >= start) & (df['start'] <= end)]
    return(gene_ref['ensembl_id'].values)

