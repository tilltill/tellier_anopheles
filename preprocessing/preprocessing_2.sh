#!/bin/bash

VCF_IN=~/vcf_files/AG1000G-CD_combined_chr1_2.vcf.gz 
VCF_OUT=~/vcf_files/AG1000G-CD_combined_chr1_2_filtered.vcf.gz 

# set filters
MAF=0.1
MISS=0.75
QUAL=30
MIN_DEPTH=20
MAX_DEPTH=100

# perform the filtering with vcftools
vcftools --gzvcf $VCF_IN \
--remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL \
--min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH \
--minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | gzip -c > \
$VCF_OUT


# Print the number of variants in the original and filtered VCF files
echo "Number of variants in the original VCF file:" 
bcftools view -H $VCF_IN | wc -l
echo "Number of variants in the filtered VCF file:"
bcftools view -H $VCF_OUT | wc -l