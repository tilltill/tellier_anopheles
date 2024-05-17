#!/usr/bin/env sh

mkdir -p $HOME/output/admixture/all_samples_QC_MAC
cd $HOME/output/admixture/all_samples_QC_MAC

# Paths
bed_file=$HOME/vcf_files/all_samples_QC_MAC.bed

# Run admixture for K = 2 to 10
for K in {2..10}; do
    admixture $bed_file $K --cv | tee log{$K}.out; done 