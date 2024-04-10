#!/bin/bash

### Download VCF files ###

# Sync the metadata using gsutil
#gsutil -m rsync -r gs://vo_agam_release/v3/metadata/ $HOME/vo_agam_release/v3/metadata/


# Set the sample_set variable
sample_sets=("AG1000G-KE")


for sample_set in "${sample_sets[@]}"; do

    # Change directory to the appropriate metadata directory
    cd ~/vo_agam_release/v3/metadata/general/$sample_set/

    # Check if the metadata CSV exists and then proceed
    csv_file="samples.meta.csv"
    if [[ -f "$csv_file" ]] && [[ $(wc -l <"$csv_file") -gt 1 ]]; then
        # Get the first five lines (header and four samples) of the metadata CSV file and save it temporarily
        head -n 26 "$csv_file" > temp_head.csv

        # Get the sample_id from the second line (first sample) of the temporary CSV file and download vcf and tbi file
        sample_id=$(awk -F, 'NR==2 {print $1}' temp_head.csv)
        wget --no-clobber -P $HOME/vcf_files "https://vo_agam_output.cog.sanger.ac.uk/$sample_id.vcf.gz" || { echo "Failed to download $sample_id.vcf.gz"; exit 1; }
        wget --no-clobber -P $HOME/vcf_files "https://vo_agam_output.cog.sanger.ac.uk/$sample_id.vcf.gz.tbi" || { echo "Failed to download $sample_id.vcf.gz"; exit 1; }
        
        # Filter first sample for chromosome 3R and 3L and save the output to a new file
        bcftools view -r 3R,3L $HOME/vcf_files/$sample_id.vcf.gz -Oz -o $HOME/vcf_files/${sample_set}_combined_chr3R3L.vcf.gz || { echo "Failed to filter $sample_id.vcf.gz"; exit 1; }

        # Create index for the new file
        bcftools index $HOME/vcf_files/${sample_set}_combined_chr3R3L.vcf.gz || { echo "Failed to index $sample_id.vcf.gz"; exit 1; }

        # Delete the original file
        rm $HOME/vcf_files/$sample_id.vcf.gz
        rm $HOME/vcf_files/$sample_id.vcf.gz.tbi

        # Read the temporary CSV file line by line starting from the third line (second sample)
        tail -n +3 temp_head.csv | while IFS=, read -r sample_id remainder; do
            # Use wget to download the .vcf.gz files for each sample_id and save them in the vcf_files directory
            wget --no-clobber -P $HOME/vcf_files "https://vo_agam_output.cog.sanger.ac.uk/$sample_id.vcf.gz" || { echo "Failed to download $sample_id.vcf.gz"; exit 1; }
            wget --no-clobber -P $HOME/vcf_files "https://vo_agam_output.cog.sanger.ac.uk/$sample_id.vcf.gz.tbi" || { echo "Failed to download $sample_id.vcf.gz"; exit 1; }

            # Create index for combined_chr3R3L.vcf.gz.csi 
            bcftools index -f $HOME/vcf_files/${sample_set}_combined_chr3R3L.vcf.gz || { echo "Failed to index ${sample_set}_combined_chr3R3L.vcf.gz"; exit 1; }

            # Merge the new file with the combined file for chromosome 3R and 3L into a temporary file
            bcftools merge $HOME/vcf_files/${sample_set}_combined_chr3R3L.vcf.gz $HOME/vcf_files/$sample_id.vcf.gz -r 3R,3L -Oz -o $HOME/vcf_files/${sample_set}_temp_combined_chr3R3L.vcf.gz || { echo "Failed to merge ${sample_set}_combined_chr3R3L.vcf.gz"; exit 1; }
            
            # Delete the original files
            rm $HOME/vcf_files/$sample_id.vcf.gz
            rm $HOME/vcf_files/$sample_id.vcf.gz.tbi

            # Move the temporary file to the original file for the next iteration
            mv $HOME/vcf_files/${sample_set}_temp_combined_chr3R3L.vcf.gz $HOME/vcf_files/${sample_set}_combined_chr3R3L.vcf.gz
        done

        rm $HOME/vcf_files/${sample_set}_temp_combined_chr3R3L.vcf.gz
        rm $HOME/vcf_files/${sample_set}_temp_combined_chr3R3L.vcf.gz.csi
    else
        echo "The metadata CSV file does not exist or is empty."
        exit 1
    fi

done