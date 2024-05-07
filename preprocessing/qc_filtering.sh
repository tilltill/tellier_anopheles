#!/bin/bash

# Path to in and out VCF files
VCF_IN=~/vcf_files/all_samples_PASS_trimmed_biallelic.vcf.gz
VCF_OUT=~/vcf_files/all_samples_QC_noMAC.vcf.gz

# Filter parameters
MAC=1
MISS=0.8
MIN_DEPTH=10
MAX_DEPTH=60

# Filter the VCF file with MAC
vcftools --gzvcf $VCF_IN \
    --mac $MAC --max-missing $MISS \
    --min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH \
    --minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | gzip -c > \
    $VCF_OUT

# Filter the VCF file without MAC
vcftools --gzvcf $VCF_IN \
    --max-missing $MISS --min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH \
    --minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | gzip -c > \
    $VCF_OUT