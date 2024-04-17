#!/bin/bash

# Path to the VCF file
input_file=~/vcf_files/AG1000G-CD_combined_chr1_2.vcf.gz 
filename=$(basename -- "$input_file")
filename="${filename%.*}"
filename="${filename%.*}"

# Make directory for the subset
mkdir -p ~/vcf_files/${filename}_analysis
mkdir -p ~/vcf_files/${filename}_analysis/vcftools

# Change to the directory
cd ~/vcf_files/${filename}_analysis

# Create random subsample of the VCF file
bcftools view $input_file | vcfrandomsample -r 0.012 > ${filename}_subset.vcf

# Compress and index subset
bcftools view ${filename}_subset.vcf -Oz -o ${filename}_subset.vcf.gz
bcftools index ${filename}_subset.vcf.gz

# Remove indels and decompose the variants
bcftools norm -m -any | \
bcftools view -V indels ${filename}_subset.vcf.gz | \
bcftools view -m2 -M2 -v snps -Oz -o ${filename}_subset.vcf.gz

# Index the new file (no indels and decomposed)
bcftools index -f ${filename}_subset.vcf.gz

# Name of the subset VCF file
SUBSET_VCF=~/vcf_files/${filename}_analysis/${filename}_subset.vcf.gz
OUT=~/vcf_files/${filename}_analysis/vcftools/${filename}_subset

# Calculate allele frequencies
vcftools --gzvcf $SUBSET_VCF --freq2 --out $OUT --max-alleles 2

# Calculate mean depth per individual
vcftools --gzvcf $SUBSET_VCF --depth --out $OUT

# Calculate mean depth per site
vcftools --gzvcf $SUBSET_VCF --site-mean-depth --out $OUT

# Calculate site quality
vcftools --gzvcf $SUBSET_VCF --site-quality --out $OUT

# Calculate missingness per individual
vcftools --gzvcf $SUBSET_VCF --missing-indv --out $OUT

# Calculate missingness per site
vcftools --gzvcf $SUBSET_VCF --missing-site --out $OUT

# Calculate heterozygosity and inbreeding coefficient per individual
vcftools --gzvcf $SUBSET_VCF --het --out $OUT

