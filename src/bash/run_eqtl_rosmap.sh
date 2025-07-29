#! /bin/bash

conda activate scanp 
cells=("Astrocyte" "Endothelial" "Microglia" "Excitatory Neuron" "Oligodendrocyte" "OPC" "Peripheral Blood" "Inhibitory Neuron" "Other")

snps=("chr2_127033851_127233851" "chr19_44844779_45044779")

for cell in "${cells[@]}"; do

    python src/run_eqtl.py -s data/adsp_r4_rosmap/chr2.rosmap.r4.wgs.biallelic.genotypes.mac10.corrected.txt.gz \
            -a data/rnaseq/rosmap_raw_rnaseq_labeled.h5ad \
            -v data/metadata/ROSMAP_covariates.csv \
            -k subject_key \
            --celltype-column major_cell_type \
            -c $cell \
            -o output/${cell}_chr2_127033851_127233851_rosmap_adsp.csv \
            --snp-range chr2_127033851_127233851
done