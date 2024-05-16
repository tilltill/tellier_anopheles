#!/bin/bash

#input_files=("$HOME/vcf_files/all_samples_QC_MAC.vcf.gz" "$HOME/vcf_files/all_samples_QC_noMAC.vcf.gz")
input_files=("$HOME/vcf_files/all_samples_QC_MAC.vcf.gz")

# Change to the directory
cd ~/vcf_files

for input_file in ${input_files[@]}; do
    # Path to the VCF file
    filename=$(basename "$input_file" .vcf.gz)

    # Edit ID field in the VCF file (otherwise IDs in prune.in will be "." which does not allow extraction after pruning)
    bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' ${filename}.vcf.gz -Oz -o ${filename}_annotated.vcf.gz

    # Convert the VCF file to PLINK format
    plink --vcf ${filename}_annotated.vcf.gz --make-bed --out ${filename}

    # Remove SNPs in LD with r^2 > 0.1 in a sliding window of 50 SNPs and step of 10 SNPs
    plink --bfile ${filename} --indep-pairwise 50 10 0.1 

    # Extract the pruned SNPs
    plink --bfile ${filename} --extract plink.prune.in --make-bed --out ${filename}_prunedData

    # Convert the pruned data back to vcf.gz format
    plink --bfile ${filename}_prunedData --recode vcf-iid bgz --out ${filename}_prunedData

    # Index the pruned VCF file
    bcftools index ${filename}_prunedData.vcf.gz

    # Print the number of variants in the original and filtered VCF files
    echo "Number of variants in the original VCF file:"
    bcftools view -H ${filename}.vcf.gz | wc -l
    echo "Number of variants in the pruned VCF file:"
    bcftools view -H ${filename}_prunedData.vcf.gz | wc -l


done

