sample_set_names=("AG1000G-CD", "AG1000G-GM-A", "AG1000G-GM-B", "AG1000G-GM-C", "AG1000G-KE",
                  "AG1000G-TZ", "1264-VO-CD-WATSENGA-VMF00161", "1264-VO-CD-WATSENGA-VMF00164", "1272-VO-GM-OPONDO-VMF00160",
                  "1246-VO-TZ-KABULA-VMF00197")

for sample_set in "${sample_set_names[@]}"; do
    
    # Create the directory structure for metadata for each sample set
    mkdir -pv $HOME/vo_agam_release/v3/metadata/general/$sample_set/

    # Download the metadata for each sample set from https://storage.googleapis.com/vo_agam_release/v3/metadata/general/$sample_set/samples.meta.csv
    wget -P $HOME/vo_agam_release/v3/metadata/general/$sample_set/ "https://storage.googleapis.com/vo_agam_release/v3/metadata/general/$sample_set/samples.meta.csv"

done
