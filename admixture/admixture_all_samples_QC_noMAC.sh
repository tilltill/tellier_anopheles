#!/usr/bin/env sh

mkdir -p $HOME/output/admixture/all_samples_QC_noMAC
cd $HOME/output/admixture/all_samples_QC_noMAC

# Paths
bed_file=$HOME/vcf_files/all_samples_QC_noMAC.bed

# Run admixture for K = 2 to 7
for K in {2..7}; do
    admixture $bed_file $K --cv | tee log{$K}.out; done 