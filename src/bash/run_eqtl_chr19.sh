#! /bin/bash

conda activate scanp 

snps=("chr2_127033851_127233851" "chr19_44844779_45044779")

cell=$1

python src/run_eqtl.py -s data/adsp_r4_rosmap/chr19.rosmap.r4.wgs.biallelic.genotypes.mac10.corrected.txt.gz \
        -a data/rnaseq/rosmap_raw_rnaseq_labeled.h5ad \
        -v data/metadata/ROSMAP_covariates.csv \
        -k subject_key \
        --celltype-column major_cell_type \
        -c "$cell" \
        -o "/project/guomiclab/hamdanz_projects/adsp_r4_rosmap/output/parallel_job_runs/${cell}_chr19_44844779_45044779_rosmap_adsp.csv" \
        --snp-range chr19_44844779_45044779 \
        --gene-range chr19_44444779_45444779


