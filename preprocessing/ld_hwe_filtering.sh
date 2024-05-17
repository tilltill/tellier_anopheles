#!/bin/bash

input_files=("$HOME/vcf_files/all_samples_QC_MAC.vcf.gz" "$HOME/vcf_files/all_samples_QC_noMAC.vcf.gz")


# Change to the directory
cd ~/vcf_files

# LD Filtering
for input_file in ${input_files[@]}; do
    # Path to the VCF file
    filename=$(basename "$input_file" .vcf.gz)

    # Edit ID field in the VCF file (otherwise IDs in prune.in will be "." which does not allow extraction after pruning)
    bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' ${filename}.vcf.gz -Oz -o ${filename}_annotated.vcf.gz

    # Convert the VCF file to PLINK format
    plink --vcf ${filename}_annotated.vcf.gz --make-bed --out ${filename}

    # LD filtering
    # Remove SNPs in LD with r^2 > 0.1 in a sliding window of 50 SNPs and step of 10 SNPs
    plink --bfile ${filename} --indep-pairwise 50 10 0.1 --out ${filename}

    # Extract the pruned SNPs
    plink --bfile ${filename} --extract ${filename}.prune.in --make-bed --out ${filename}_LD

    # HWE filtering for SNPs with HWE p-value < 1e-50 for LD pruned data
    plink --bfile ${filename}_LD --hwe 1e-50 --make-bed --out ${filename}_LD_HWE

    # HWE filtering for SNPs with HWE p-value < 1e-50 for non-LD pruned data
    plink --bfile ${filename} --hwe 1e-50 --make-bed --out ${filename}_HWE

    # Convert filtered data back to vcf.gz format
    plink --bfile ${filename}_LD --recode vcf-iid bgz --out ${filename}_LD
    plink --bfile ${filename}_LD_HWE --recode vcf-iid bgz --out ${filename}_LD_HWE
    plink --bfile ${filename}_HWE --recode vcf-iid bgz --out ${filename}_HWE

    # Index the filtered VCF files
    bcftools index ${filename}_LD.vcf.gz
    bcftools index ${filename}_LD_HWE.vcf.gz
    bcftools index ${filename}_HWE.vcf.gz


    # Print the number of variants in the original LD pruned, HWE pruned and LD+HWE pruned VCF files
    echo "Number of variants in the original VCF file:"
    bcftools view -H ${filename}.vcf.gz | wc -l
    
    echo "Number of variants in the LD pruned VCF file:"
    bcftools view -H ${filename}_LD.vcf.gz | wc -l
    
    echo "Number of variants in the HWE pruned VCF file:"
    bcftools view -H ${filename}_HWE.vcf.gz | wc -l
    
    echo "Number of variants in the LD+HWE pruned VCF file:"
    bcftools view -H ${filename}_LD_HWE.vcf.gz | wc -l
    
done