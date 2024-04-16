#!/bin/bash

input_file=

# Get filename without extension or path
filename=$(basename -- "$input_file")
filename="${filename%.*}"

# Remove indels 
bcftools view -V indels $input_file -Oz -o ${filename}_preprocessed.vcf.gz

# Index the new file
bcftools index -f ${filename}_preprocessed.vcf.gz

# Decompose the variants
bcftools norm -m -any ${filename}_preprocessed.vcf.gz -Oz -o ${filename}_preprocessed.vcf.gz

# Index the new file
bcftools index -f ${filename}_preprocessed.vcf.gz

# Filter for biallelic SNPs
bcftools view -m2 -M2 -v snps ${filename}_preprocessed.vcf.gz -Oz -o ${filename}_preprocessed.vcf.gz

# Index the new file
bcftools index -f ${filename}_preprocessed.vcf.gz

# Add missingness Info
bcftools +fill-tags $input_file -Oz -o ${filename}_preprocessed.vcf.gz -- -t F_MISSING

# Index the new file
bcftools index -f ${filename}_preprocessed.vcf.gz

# Remove variants with F_MISSING > 0.25
bcftools view ${filename}_preprocessed.vcf.gz -i 'F_MISSING < 0.3' -Oz -o ${filename}_preprocessed.vcf.gz

# Index the new file
bcftools index -f ${filename}_preprocessed.vcf.gz

# Remove variants with read depth < 10 or > 100
bcftools filter -i 'AVG(DP)>10 && AVG(DP)<100' ${filename}_preprocessed.vcf.gz -Oz -o ${filename}_preprocessed.vcf.gz

# Index the new file
bcftools index -f ${filename}_preprocessed.vcf.gz

# Remove variants with Phred score < 30
bcftools filter -i 'AVG(GQ)>30' ${filename}_preprocessed.vcf.gz -Oz -o ${filename}_preprocessed.vcf.gz

# Index the new file
bcftools index -f ${filename}_preprocessed.vcf.gz

