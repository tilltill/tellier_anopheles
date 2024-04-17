#!/bin/bash

# Path to the VCF file
input_file=~/vcf_files/AG1000G-CD_combined_chr1_2.vcf.gz 

# Change to the directory
cd ~/vcf_files

# Get filename without extension or path
filename=$(basename -- "$input_file")
filename="${filename%.*}"

# Chain commands using pipes
bcftools view -V indels $input_file | \
bcftools norm -m -any | \
bcftools view -m2 -M2 -v snps | \
bcftools +fill-tags -- -t F_MISSING | \
bcftools view -i 'F_MISSING < 0.25' | \
bcftools filter -i 'AVG(DP)>10 && AVG(DP)<100' | \
bcftools filter -i 'AVG(GQ)>30' -Oz -o ${filename}_preprocessed.vcf.gz

# Index the final file
bcftools index -f ${filename}_preprocessed.vcf.gz