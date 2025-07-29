#! /bin/bash

conda activate scanp 

snps=("chr2_127033851_127233851" "chr19_44844779_45044779")

cell=$1

python src/run_eqtl.py -s /project/guomiclab/hamdanz_projects/adsp_r4_rosmap/data/adsp_r4_rosmap/chr2.rosmap.r4.wgs.biallelic.genotypes.mac10.corrected.txt.gz \
        -a data/rnaseq/rosmap_raw_rnaseq_labeled.h5ad \
        -v data/metadata/ROSMAP_covariates.csv \
        -k subject_key \
        --celltype-column major_cell_type \
        -c "$cell" \
        -o "/project/guomiclab/hamdanz_projects/adsp_r4_rosmap/eqtl_results_05222025/chr2_127033851_127233851_gene_range_126633851_127633851/${cell}_chr2_127033851_127233851_rosmap_adsp.csv" \
        --snp-range chr2_127033851_127233851 \
        --gene-range chr2_126633851_127633851 