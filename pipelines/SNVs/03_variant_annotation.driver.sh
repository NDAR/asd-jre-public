#!/bin/bash

# test if the first argument is unset, if not, set it to the FAMILY_ID variable
if [ ! -z "$1" ]; 
then
    FAMILY_ID=$1
    echo "FAMILY_ID is ${FAMILY_ID}"
fi


SCRIPT_DIR="03_variant_annotation_pipeline"
MANIFEST_PATH="/path/to/manifest.txt"

PROJECT_ID="default_project"
HOSTNAME=`hostname`

# S3 bucket name for storing data
S3_BUCKET="s3://your-bucket-name"

echo "Starting annotation_driver.sh for ${FAMILY_ID}"

mkdir -p /mnt/reference
mkdir -p /mnt/temp
mkdir -p /mnt/out

s3cmd get --skip-existing $S3_BUCKET/REFERENCE/annotation/* /mnt/reference/
s3cmd get --skip-existing $S3_BUCKET/REFERENCE/dbsnp_137.b37.vcf.gz /mnt/reference/dbsnp_137.b37.vcf

# note we do not gunzip everything because compresed vcfs are used by annotation tools
# TODO check if these are unzipped already
if [ ! -e /mnt/reference/dbNSFP2.1.txt ]; then
    gunzip -c /mnt/reference/dbNSFP2.1.txt.gz > /mnt/reference/dbNSFP2.1.txt
fi
if [ ! -e /mnt/reference/cadd_v1.0.nimblegen_solution_V2refseq2010.HG19.300bp_slop.vcf ]; then
    gunzip -c /mnt/reference/cadd_v1.0.nimblegen_solution_V2refseq2010.HG19.300bp_slop.vcf.gz > /mnt/reference/cadd_v1.0.nimblegen_solution_V2refseq2010.HG19.300bp_slop.vcf
fi
if [ ! -e /mnt/reference/dbsnp_137.b37.vcf ]; then
    gunzip -c /mnt/reference/dbsnp_137.b37.vcf.gz > /mnt/reference/dbsnp_137.b37.vcf
fi

apt-get install --force-yes -q -y tabix
export PERL5LIB=/usr/bin/vcftools/vcftools_0.1.11/perl

#######
# VariantAnnotation_Pipeline
# 
python ${SCRIPT_DIR}/VariantAnnotation_Pipeline.py ${RUFFUS_ARGS} ${MANIFEST_PATH} ${FAMILY_ID} /mnt ${OPTS}
if [ $? -ne 0 ]; then
    echo "Non-zero return code received from VariantAnnotation_Pipeline.py"
    exit 1
fi

echo "Completed VariantAnnotation_Pipeline for ${FAMILY_ID}"
