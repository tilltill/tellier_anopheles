#!/bin/bash


# Set working directory
cd ~/vcf_files

# Kenyan samples
echo STARTED merging Kenyan samples on $(date)
bcftools merge AG1000G-KE_combined_chr1_2.vcf.gz chr1_sitefilter.vcf.gz -r 1 -Oz -0 AG1000G-KE_combined_chr1_2_sitefilter_chr1.vcf.gz
bcftools index -f AG1000G-KE_combined_chr1_2_sitefilter_chr1.vcf.gz

bcftools merge AG1000G-KE_combined_chr1_2_sitefilter_chr1.vcf.gz chr2_sitefilter.vcf.gz -r 2 -Oz -o AG1000G-KE_combined_chr1_2_sitefilter_chr1_2.vcf.gz
bcftools index -f AG1000G-KE_combined_chr1_2_sitefilter_chr1_2.vcf.gz

rm -f AG1000G-KE_combined_chr1_2_sitefilter_chr1.vcf.gz AG1000G-KE_combined_chr1_2_sitefilter_chr1.vcf.gz.csi

echo FINISHED merging Kenyan samples on $(date)

# Gambian samples
echo STARTED merging Gambian samples on $(date)
bcftools merge AG1000G-GM-C_combined_chr1_2.vcf.gz chr1_sitefilter.vcf.gz -r 1 -Oz -o AG1000G-GM-C_combined_chr1_2_sitefilter_chr1.vcf.gz
bcftools index -f AG1000G-GM-C_combined_chr1_2_sitefilter_chr1.vcf.gz

bcftools merge AG1000G-GM-C_combined_chr1_2_sitefilter_chr1.vcf.gz chr2_sitefilter.vcf.gz -r 2 -Oz -o AG1000G-GM-C_combined_chr1_2_sitefilter_chr1_2.vcf.gz
bcftools index -f AG1000G-GM-C_combined_chr1_2_sitefilter_chr1_2.vcf.gz

rm -f AG1000G-GM-C_combined_chr1_2_sitefilter_chr1.vcf.gz AG1000G-GM-C_combined_chr1_2_sitefilter_chr1.vcf.gz.csi

echo FINISHED merging Gambian samples on $(date)

# Tanzanian samples
echo STARTED merging Tanzanian samples on $(date)
bcftools merge AG1000G-TZ_combined_chr3R3L.vcf.gz chr1_sitefilter.vcf.gz -r 1 -Oz -o AG1000G-TZ_combined_chr3R3L_sitefilter_chr1.vcf.gz
bcftools index -f AG1000G-TZ_combined_chr3R3L_sitefilter_chr1.vcf.gz

bcftools merge AG1000G-TZ_combined_chr3R3L_sitefilter_chr1.vcf.gz chr2_sitefilter.vcf.gz -r 2 -Oz -o AG1000G-TZ_combined_chr3R3L_sitefilter_chr1_2.vcf.gz
bcftools index -f AG1000G-TZ_combined_chr3R3L_sitefilter_chr1_2.vcf.gz

rm -f AG1000G-TZ_combined_chr3R3L_sitefilter_chr1.vcf.gz AG1000G-TZ_combined_chr3R3L_sitefilter_chr1.vcf.gz.csi

echo FINISHED merging Tanzanian samples on $(date)

