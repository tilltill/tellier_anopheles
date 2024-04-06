#!/bin/bash

### Download VCF files ###

# Set the sample_set variable 
sample_set="AG1000G-CD"

# Change directory to the appropriate metadata directory
cd ~/vo_agam_release/v3/metadata/general/$sample_set/

# Check if the metadata CSV exists and then proceed
csv_file="samples.meta.csv"
if [[ -f "$csv_file" ]] && [[ $(wc -l <"$csv_file") -gt 1 ]]; then
    # Get the first five lines (header and four samples) of the metadata CSV file and save it temporarily
    head -n 5 "$csv_file" > temp_head.csv

    # Get the sample_id from the second line (first sample) of the temporary CSV file
    sample_id=$(awk -F, 'NR==2 {print $1}' temp_head.csv)
    wget --no-clobber -P $HOME/vcf_files/test_download "https://vo_agam_output.cog.sanger.ac.uk/$sample_id.vcf.gz" || { echo "Failed to download $sample_id.vcf.gz"; exit 1; }
    wget --no-clobber -P $HOME/vcf_files/test_download "https://vo_agam_output.cog.sanger.ac.uk/$sample_id.vcf.gz.tbi" || { echo "Failed to download $sample_id.vcf.gz"; exit 1; }

    # Filter first sample for chromosome 3R and save the output to a new file
    bcftools view -r 3R $HOME/vcf_files/test_download/$sample_id.vcf.gz -Oz -o $HOME/vcf_files/test_download/combined_chr3R.vcf.gz || { echo "Failed to filter $sample_id.vcf.gz"; exit 1; }

    # Create index for the new file
    bcftools index $HOME/vcf_files/test_download/combined_chr3R.vcf.gz || { echo "Failed to index $sample_id.vcf.gz"; exit 1; }

    # Delete the original file
    #rm $HOME/vcf_files/test_download/$sample_id.vcf.gz
    #rm $HOME/vcf_files/test_download/$sample_id.vcf.gz.tbi

    # Read the temporary CSV file line by line starting from the third line (second sample)
    tail -n +3 temp_head.csv | while IFS=, read -r sample_id remainder; do
        # Use wget to download the .vcf.gz files for each sample_id and save them in the vcf_files directory
        wget --no-clobber -P $HOME/vcf_files/test_download "https://vo_agam_output.cog.sanger.ac.uk/$sample_id.vcf.gz" || { echo "Failed to download $sample_id.vcf.gz"; exit 1; }
        wget --no-clobber -P $HOME/vcf_files/test_download "https://vo_agam_output.cog.sanger.ac.uk/$sample_id.vcf.gz.tbi" || { echo "Failed to download $sample_id.vcf.gz"; exit 1; }

        # Filter for chromosome 3R and save the output to a new file
        bcftools view -r 3R $HOME/vcf_files/test_download/$sample_id.vcf.gz -Oz -o $HOME/vcf_files/test_download/${sample_id}_chr3R.vcf.gz || { echo "Failed to filter $sample_id.vcf.gz"; exit 1; }
        
        # Create index for the new filtered file
        bcftools index $HOME/vcf_files/test_download/${sample_id}_chr3R.vcf.gz || { echo "Failed to index ${sample_id}_chr3R.vcf.gz"; exit 1; }

        # Merge the new file with the combined file for chromosome 3R into a temporary file
        bcftools merge $HOME/vcf_files/test_download/combined_chr3R.vcf.gz $HOME/vcf_files/test_download/${sample_id}_chr3R.vcf.gz  -Oz -o $HOME/vcf_files/test_download/temp_combined_chr3R.vcf.gz || { echo "Failed to merge combined_chr3R.vcf.gz and ${sample_id}_chr3R.vcf.gz"; exit 1; }

        # Create index for the temp file
        bcftools index -f $HOME/vcf_files/test_download/temp_combined_chr3R.vcf.gz -o $HOME/vcf_files/test_download/combined_chr3R.vcf.gz.csi || { echo "Failed to index $sample_id.vcf.gz"; exit 1; }

        # Delete the original files
        #rm $HOME/vcf_files/test_download/$sample_id.vcf.gz
        #rm $HOME/vcf_files/test_download/$sample_id.vcf.gz.tbi
        #rm $HOME/vcf_files/test_download/${sample_id}_chr3R.vcf.gz

        # Move the temporary file to the original file for the next iteration
        mv $HOME/vcf_files/test_download/temp_combined_chr3R.vcf.gz $HOME/vcf_files/test_download/combined_chr3R.vcf.gz
    done

    #rm temp_head.csv
    #rm $HOME/vcf_files/test_download/temp_combined_chr3R.vcf.gz
    #rm $HOME/vcf_files/test_download/temp_combined_chr3R.vcf.gz.csi
else
    echo "The metadata CSV file does not exist or is empty."
    exit 1
fi