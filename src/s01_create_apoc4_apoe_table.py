import scanpy as sc
import pandas as pd

adata_path = '/project/guomiclab/hamdanz_projects/adsp_r4_rosmap/data/rnaseq/rosmap_raw_rnaseq_labeled.h5ad'
snp_path = '/project/guomiclab/hamdanz_projects/adsp_r4_rosmap/data/adsp_r4_rosmap/chr19.rosmap.r4.wgs.biallelic.genotypes.mac10.txt.gz'

target_snp = 'chr19_44944779_T_C'

adata = sc.read_h5ad(adata_path)

snp_df = pd.read_csv(snp_path, sep='\t')
snp_df.columns = [col.split('-B')[0] for col in snp_df.columns]
snp_df.columns = [col.split('-U')[0] for col in snp_df.columns]
snp_df['snp_id'] = snp_df['CHROM'].astype(str) + '_' + snp_df['POS'].astype(str) + '_' + snp_df['REF'].astype(str) + '_' + snp_df['ALT'].astype(str)
snp_df = snp_df[snp_df['snp_id'] == target_snp]
snp_df.index = snp_df['snp_id']

adata.var.index = adata.var.gene_name
adata = adata[:, ~adata.var.index.duplicated(keep='last')]
adata = adata[adata.obs.major_cell_type == 'Microglia']

overlap_ids = list(set(adata.obs.subject_key).intersection(set(snp_df.columns)))
snp_df = snp_df[overlap_ids]
adata= adata[adata.obs.subject_key.isin(overlap_ids)]

apoe_apoc_adata = adata[:,['APOC4', 'APOE']]
mglia_apoe_apoc_adata = apoe_apoc_adata[apoe_apoc_adata.obs.major_cell_type == 'Microglia']
agg_apo_adata = sc.get.aggregate(mglia_apoe_apoc_adata, by='individualID', func='mean')

individual_snp_id_match = mglia_apoe_apoc_adata.obs[['subject_key', 'individualID']]
individual_snp_id_match = individual_snp_id_match.drop_duplicates()
individual_snp_id_match.index = individual_snp_id_match.individualID

mglia_agg_expr_df = pd.DataFrame(agg_apo_adata.layers['mean'])
mglia_agg_expr_df.columns = list(agg_apo_adata.var.index.astype(str))
mglia_agg_expr_df.index = agg_apo_adata.obs.index.astype(str)
mglia_agg_expr_df = pd.merge(mglia_agg_expr_df, individual_snp_id_match, left_index=True, right_index=True)
mglia_agg_expr_df.index = mglia_agg_expr_df['subject_key'].astype(str)
mglia_agg_expr_df = mglia_agg_expr_df.drop(columns=['subject_key','individualID']) 
log_mglia_agg_expr_df = mglia_agg_expr_df.copy()
log_mglia_agg_expr_df = log_mglia_agg_expr_df.apply(lambda x: np.log(x + 1))

snp_genotype_match = pd.merge(mglia_agg_expr_df, snp_df.T, left_index=True, right_index=True)
snp_genotype_match.index = snp_genotype_match.index.astype(str)

lognorm_snp_genotype_match = pd.merge(log_mglia_agg_expr_df, snp_df.T, left_index=True, right_index=True)
lognorm_snp_genotype_match.index = snp_genotype_match.index.astype(str)

snp_genotype_match= snp_genotype_match[['chr19_44944779_T_C', 'APOC4', 'APOE']]
lognorm_snp_genotype_match = lognorm_snp_genotype_match[['chr19_44944779_T_C', 'APOC4', 'APOE']]

snp_genotype_match.to_csv('/project/guomiclab/hamdanz_projects/adsp_r4_rosmap/output/apoc4_apoe_chr19_44944779_T_C_mean_agg_expr_mglia_only.csv', index=True)
lognorm_snp_genotype_match.to_csv('/project/guomiclab/hamdanz_projects/adsp_r4_rosmap/output/apoc4_apoe_chr19_44944779_T_C_lognorm_mean_agg_expr_mglia_only.csv', index=True)