#!/usr/bin/env sh

cd $HOME/output/admixture/all_samples

# Paths
bed_file=$HOME/vcf_files/all_samples_PASS_trimmed_biallelic.bed

# Run admixture for K = 2 to 10
for K in {2..10}; do
    admixture $bed_file $K --cv | tee log{$K}.out; done 