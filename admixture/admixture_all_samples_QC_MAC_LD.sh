#!/usr/bin/env sh

mkdir -p $HOME/output/admixture/all_samples_QC_MAC_LD
cd $HOME/output/admixture/all_samples_QC_MAC_LD

# Paths
bed_file=$HOME/vcf_files/all_samples_QC_MAC_LD.bed

# Run admixture for K = 2 to 7
for K in {2..7}; do
    admixture $bed_file $K --cv | tee log{$K}.out; done 