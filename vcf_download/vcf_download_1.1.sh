#!/bin/bash

#SBATCH --job-name=vcf_anopheles_download_test1
#SBATCH --output=vcf_anopheles_download_test1.out
#SBATCH --error=vcf_anopheles_download_test1.err
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G


### Download VCF files ###

# Set the sample_set variable
sample_set="AG1000G-CD"

# gsutil not installed on server, so I uploaded the metadata files to the server manually

# Create the directory structure for metadata
#mkdir -pv $HOME/vo_agam_release/v3/metadata/

# Create the directory structure for VCF files
#mkdir -p $HOME/vcf_files


# Sync the metadata using gsutil
#gsutil -m rsync -r gs://vo_agam_release/v3/metadata/ $HOME/vo_agam_release/v3/metadata/

# Change directory to the appropriate metadata directory
cd ~/vo_agam_release/v3/metadata/general/$sample_set/

# Check if the metadata CSV exists and if it has more than one line (excluding the header) then proceed
csv_file="samples.meta.csv"
if [[ -f "$csv_file" ]] && [[ $(wc -l <"$csv_file") -gt 1 ]]; then
    # Get the first eleven lines of the metadata CSV file and save it temporarily
    head -n 11 "$csv_file" > temp_head.csv

    # Initialize an empty VCF file for concatenation
    touch $HOME/vcf_files/combined_chr3.vcf.gz

    # Read the temporary CSV file line by line, excluding the header
    tail -n +2 temp_head.csv | while IFS=, read -r sample_id remainder; do
    
        # Use wget to download the .vcf.gz files for each sample_id and save them in the vcf_files directory
        wget --no-clobber -P $HOME/vcf_files "https://vo_agam_output.cog.sanger.ac.uk/$sample_id.vcf.gz" || { echo "Failed to download $sample_id.vcf.gz"; exit 1; }

        # Download the corresponding index file
        #wget --no-clobber -P $HOME/vcf_files "https://vo_agam_output.cog.sanger.ac.uk/$sample_id.vcf.gz.tbi" || { echo "Failed to download $sample_id.vcf.gz.tbi"; exit 1; }
        
        # Merge the new file with the combined file for chromosome 3
        bcftools merge $HOME/vcf_files/combined_chr3.vcf.gz $HOME/vcf_files/$sample_id.vcf.gz -Oz -r 3 -o $HOME/vcf_files/combined_chr3.vcf.gz
        
        # Delete the original file
        rm $HOME/vcf_files/$sample_id.vcf.gz

        # Concatenate the new file with the combined file
        #bcftools concat $HOME/vcf_files/combined_chr3.vcf.gz $HOME/vcf_files/${sample_id}_chr3.vcf.gz -o $HOME/vcf_files/temp_combined.vcf.gz -O z || { echo "Failed to concatenate $sample_id.vcf.gz"; exit 1; }

        # Move the temporary file to the original file for the next iteration
        #mv $HOME/vcf_files/temp_combined.vcf.gz $HOME/vcf_files/combined_chr3.vcf.gz

        # Delete the new file
        rm $HOME/vcf_files/${sample_id}_chr3.vcf.gz
    done
    
    # Clean up the temporary file
    rm temp_head.csv
else
    echo "The metadata CSV file does not exist or is empty."
    exit 1
fi