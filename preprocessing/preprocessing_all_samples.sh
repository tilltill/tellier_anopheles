#!/bin/bash

# Path to the VCF file
input_files=("~/vcf_files/AG1000G-CD_combined_chr1_2.vcf.gz" \
             "~/vcf_files/AG1000G-GM-C_combined_chr1_2.vcf.gz" \
             "~/vcf_files/AG1000G-KE_combined_chr1_2.vcf.gz" \
             "~/vcf_files/AG1000G-TZ_combined_chr1_2.vcf.gz")

# Change to the directory
cd ~/vcf_files

for input_file in "${input_files[@]}"; do

    # Get filename without extension or path
    filename=$(basename "$input_file" .vcf.gz)

    # Merge the VCF file with the sitefilter file
    output_file=${filename}_sitefilter.vcf.gz
    if [ ! -f "$output_file" ]; then
        bcftools merge $input_file chr1and2_sitefilter.vcf.gz -Oz -o $output_file
        bcftools index -f $output_file
        echo "$(date): Merged VCF file with sitefilter file and indexed the new file"
    fi

    # Filter for PASS variants
    output_file=${filename}_sitefilter_PASS.vcf.gz
    if [ ! -f "$output_file" ]; then
        bcftools view  ${filename}_sitefilter.vcf.gz -f PASS -Oz -o $output_file
        bcftools index -f $output_file
        echo "$(date): Filtered for PASS variants and indexed the new file"
    fi

    # Trim alternative alleles
    output_file=${filename}_sitefilter_PASS_trimmed.vcf.gz
    if [ ! -f "$output_file" ]; then
        bcftools view -a ${filename}_sitefilter_PASS.vcf.gz -Oz -o $output_file
        bcftools index -f $output_file
        echo "$(date): Trimmed alternative alleles and indexed the new file"
    fi

    # Only keep biallelic SNPs
    output_file=${filename}_sitefilter_PASS_trimmed_biallelic.vcf.gz
    if [ ! -f "$output_file" ]; then
        bcftools view -m2 -M2 -v snps ${filename}_sitefilter_PASS_trimmed.vcf.gz -Oz -o $output_file
        bcftools index -f $output_file
        echo "$(date): Kept only biallelic SNPs and indexed the new file"
    fi


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

done