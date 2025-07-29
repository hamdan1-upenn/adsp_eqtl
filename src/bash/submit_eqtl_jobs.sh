#! /bin/bash

cells=("Excitatory Neuron" "Endothelial" "Oligodendrocyte" "OPC" "Peripheral Blood" "Inhibitory Neuron" "Other")

chr2_snp_file='/project/guomiclab/hamdanz_projects/adsp_r4_rosmap/data/adsp_r4_rosmap/chr2.rosmap.r4.wgs.biallelic.genotypes.mac10.corrected.txt.gz'
chr19_snp_file='/project/guomiclab/hamdanz_projects/adsp_r4_rosmap/data/adsp_r4_rosmap/chr19.rosmap.r4.wgs.biallelic.genotypes.mac10.corrected.txt.gz'

snps=("chr2_127033851_127233851" "chr19_44844779_45044779")
cells=("Excitatory Neuron")
for cell in "${cells[@]}"; do
    bsub -o "logs/${cell}_chr2_127033851_127233851_gr_126633851_127633851.out" -M 180000 "bash run_eqtl_chr2.sh \"${cell}\""
done 
for cell in "${cells[@]}"; do
    bsub -o "logs/${cell}_chr19_44844779_45044779.out" -M 180000 "bash run_eqtl_chr19.sh \"${cell}\""
done
