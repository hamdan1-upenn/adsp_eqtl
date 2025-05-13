import pandas as pd

WORK_DIR = '/project/guomiclab/hamdanz_projects/adsp_r4_rosmap'
SNP_META_PATH = f'{WORK_DIR}/data/metadata/ROSMAP_ADSPID_mapping.csv'
CLINICAL_META_PATH = f'{WORK_DIR}/data/metadata/ROSMAP_clinical.csv'
BIOSPEC_META_PATH = f'{WORK_DIR}/data/metadata/ROSMAP_biospecimen_metadata.csv'

snp_meta = pd.read_csv(SNP_META_PATH)
clinical_meta = pd.read_csv(CLINICAL_META_PATH)
biospec_meta = pd.read_csv(BIOSPEC_META_PATH)

# this string should match biospecimen metadata specimenID
snp_meta['specimenID'] = snp_meta['cohort_key'].astype(str) + snp_meta['subject_key_original'].astype(str)


merged_meta = snp_meta.merge(biospec_meta[['individualID','specimenID']], how='left', left_on='specimenID', right_on='specimenID')

# only keep individuals with corresponding individual IDs
merged_meta = merged_meta[~merged_meta.individualID.isnull()]

merged_meta = merged_meta.merge(clinical_meta[["individualID", "msex", 'age_death', "braaksc","ceradsc","cogdx","dcfdx_lv"]], left_on='individualID', right_on='individualID', how='left')

merged_meta.to_csv(f'{WORK_DIR}/data/metadata/ROSMAP_ADSPID_merged_metadata.csv', index=False)



