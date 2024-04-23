#!/bin/bash

# Sync the site filters using gsutil
mkdir -pv ~/vo_agam_release/v3/site_filters/dt_20200416/vcf/
gsutil -m rsync -r \
    gs://vo_agam_release/v3/site_filters/dt_20200416/vcf/ \
    ~/vo_agam_release/v3/site_filters/dt_20200416/vcf/

