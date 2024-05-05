#!/usr/bin/env sh

cd $HOME/output/admixture/DRC_AG1000G-CD

# Paths
bed_file=$HOME/vcf_files/AG1000G-CD_combined_chr1_2_sitefilter_PASS_trimmed_biallelic.bed

# Run admixture for K = 2 to 10
for K in {2..10}; do
    admixture $bed_file $K --cv | tee log{$K}.out; done 