#!/bin/bash

# Path to the VCF file
input_file=~/vcf_files/AG1000G-GM-C_combined_chr1_2.vcf.gz 

# Change to the directory
cd ~/vcf_files

# Get filename without extension or path
filename=$(basename "$input_file" .vcf.gz)

# Merge the VCF file with the sitefilter file
bcftools merge $input_file chr1and2_sitefilter.vcf.gz -Oz -o ${filename}_sitefilter.vcf.gz
bcftools index -f ${filename}_sitefilter.vcf.gz

# Filter for PASS variants
bcftools view  ${filename}_sitefilter.vcf.gz -f PASS -Oz -o ${filename}_sitefilter_PASS.vcf.gz
bcftools index -f ${filename}_sitefilter_PASS.vcf.gz

# Trim alternative alleles
bcftools view -a -s ${filename}_sitefilter_PASS.vcf.gz -Oz -o ${filename}_sitefilter_PASS_trimmed.vcf.gz
bcftools index -f ${filename}_sitefilter_PASS_trimmed.vcf.gz

# Only keep biallelic SNPs
bcftools view -m2 -M2 -v snps ${filename}_sitefilter_PASS_trimmed.vcf.gz -Oz -o ${filename}_sitefilter_PASS_trimmed_biallelic.vcf.gz
bcftools index -f ${filename}_sitefilter_PASS_trimmed_biallelic.vcf.gz

#bcftools +fill-tags -- -t F_MISSING | \
#bcftools view -i 'F_MISSING < 0.25' | \
#bcftools filter -i 'AVG(DP)>10 && AVG(DP)<100' | \
#bcftools filter -i 'AVG(GQ)>30' -Oz -o ${filename}_preprocessed.vcf.gz

# Print the number of variants in the original and filtered VCF files
echo "Number of variants in the original VCF file:" 
bcftools view -H $input_file | wc -l

echo "Number of variants in the filtered VCF file:"
bcftools view -H ${filename}_sitefilter_PASS.vcf.gz | wc -l

echo "Number of variants in the trimmed VCF file:"
bcftools view -H ${filename}_sitefilter_PASS_trimmed.vcf.gz | wc -l

echo "Number of variants in the biallelic VCF file:"
bcftools view -H ${filename}_sitefilter_PASS_trimmed_biallelic.vcf.gz | wc -l

