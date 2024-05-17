#!/bin/bash

input_files=("$HOME/vcf_files/all_samples_QC_MAC.vcf.gz")


# Change to the directory
cd $HOME/vcf_files

# Convert the VCF file to PLINK format and back to vcf.gz format to test if the conversion works

for input_file in ${input_files[@]}; do
    # Path to the VCF file
    filename=$(basename "$input_file" .vcf.gz)

    # Convert the VCF file to PLINK format
    plink --vcf ${filename}.vcf.gz --make-bed --out ${filename}_test

    # Convert filtered data back to vcf.gz format
    plink --bfile ${filename}_test --recode vcf-iid bgz --out ${filename}_test

    # Print the number of variants in the original LD pruned, HWE pruned and LD+HWE pruned VCF files
    echo "Number of variants in the original VCF file:"
    bcftools view -H ${filename}.vcf.gz | wc -l
    
    echo "Number of variants in the test file:"
    bcftools view -H ${filename}_test.vcf.gz | wc -l

done