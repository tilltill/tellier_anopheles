#!/bin/bash


### Download VCF files 
# Set the sample_set variable
sample_set="AG1000G-CD"

# Create the directory structure for metadata
mkdir -pv ~/vo_agam_release/v3/metadata/

# Sync the metadata using gsutil
gsutil -m rsync -r gs://vo_agam_release/v3/metadata/ ~/vo_agam_release/v3/metadata/

# Change directory to the appropriate metadata directory
cd ~/vo_agam_release/v3/metadata/general/$sample_set/

# Check if the metadata CSV exists and then proceed
csv_file="samples.meta.csv"
if [[ -f "$csv_file" ]]; then
    # Get the first ten lines of the metadata CSV file and save it temporarily
    head -n 10 "$csv_file" > temp_head.csv
    
    # Create output directory for vcf files (if it does not exist)
    mkdir -p $HOME/vcf_files

    # Read the temporary CSV file line by line, excluding the header
    tail -n +2 temp_head.csv | while IFS=, read -r sample_id remainder; do
        # Use wget to download the .vcf.gz files for each sample_id and save them in the vcf_files directory
        wget --no-clobber -P $HOME/vcf_files "https://vo_agam_output.cog.sanger.ac.uk/$sample_id.vcf.gz"
    done
    
    # Clean up the temporary file
    rm temp_head.csv
else
    echo "The metadata CSV file does not exist."
    exit 1
fi

### Combine VCF files for chromosome x ###
cd $HOME/vcf_files

# Step 1: Filter for chromosome 1
for file in *.vcf; do
    bcftools view -r 1 $file -o ${file%.vcf}_chr1.vcf
done

# Step 2: Combine filtered files
bcftools concat *_chr1.vcf -o combined_chr1.vcf -O v

# Step 3: Delete original files (be careful!)
rm *.vcf
