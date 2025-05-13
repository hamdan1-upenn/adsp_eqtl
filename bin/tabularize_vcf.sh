#! /bin/bash

bcftools_path=$1
vcf_path=$2

tmp_dir=$(dirname $vcf_path)/tmp

mkdir -p $tmp_dir

if [ -z "$bcftools_path" ] || [ -z "$vcf_path" ]; then
    echo "Usage: $0 <bcftools_path> <vcf_path>"
    exit 1
fi

(echo -e "CHROM\nPOS\nREF\nALT\nRSID\nAC\nAN\nAF"; $bcftools_path query -l $vcf_path) | tr '\n' '\t'; 
    $bcftools_path view -c 10 -i 'INFO/AN > 1314' $vcf_path  | \
    $bcftools_path query -f '\n%CHROM\t%POS\t%REF\t%ALT\t%ID\t%INFO/AC\t%INFO/AN\t%INFO/AF[\t%GT]'   | \
    sed -e 's/0\/0/0/g' -e 's/0\/1/1/g' -e 's/1\/0/1/g' -e 's/1\/1/2/g' -e 's/.\/./NA/g' -e 's/0|0/0/g' -e 's/0|1/1/g' -e 's/1|0/1/g' -e 's/1|1/2/g' -e 's/.\/./NA/g' \
    > $tmp_dir/tabularized_vcf.tsv
