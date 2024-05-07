#!/bin/bash

# Path to the VCF file and filename
input_file=~/vcf_files/all_samples_PASS_trimmed_biallelic.vcf.gz
filename=$(basename "$input_file" .vcf.gz)

# Filter parameters
MAC=1
MISS=0.8
MIN_DEPTH=10
MAX_DEPTH=60

# Filter the VCF file with MAC
vcftools --gzvcf ${filename}_sitefilter_PASS_trimmed.vcf.gz \
    --mac $MAC --max-missing $MISS \
    --min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH \
    --minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | gzip -c > \
    ${filename}_QC_MAC.vcf.gz

# Filter the VCF file without MAC
vcftools --gzvcf ${filename}_sitefilter_PASS_trimmed.vcf.gz \
    --max-missing $MISS --min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH \
    --minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | gzip -c > \
    ${filename}_QC_noMAC.vcf.gz