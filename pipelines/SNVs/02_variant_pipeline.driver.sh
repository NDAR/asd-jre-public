#!/bin/bash

# test if the first argument is unset, if not, set it to the BATCH_ID variable
if [ ! -z "$1" ]; 
then
    BATCH_ID=$1
    echo "BATCH_ID is ${BATCH_ID}"
fi

SCRIPT_DIR="02_variant_pipeline"
MANIFEST_PATH="/path/to/batches.txt"

PROJECT_ID="default_project"
HOSTNAME=`hostname`

FAMILY_ID_LIST=$(grep ${BATCH_ID} ${MANIFEST_PATH} | cut -f2 -d, | sort | uniq)

echo "Starting variant_driver (batch ${BATCH_ID})"

####### 
# mkdirs as appropriate
#
mkdir -p /mnt/temp
mkdir -p /mnt/reference
mkdir -p /mnt/out

#######
# - Download reference files
# - gunzip them
s3cmd get --recursive --skip-existing s3://asdjre/REFERENCE/ /mnt/reference/
gunzip -f /mnt/reference/*.gz


#######
# Run GATK_variants_batch.py script
# - downloads reduced read bams from S3
# - Runs HC
# - Calcualtes and applies VQSR
# - splits up VCF file and uploads to each output dir on S3

if [ ! -n "$SKIP_GATK" ]; then

    python ${SCRIPT_DIR}/GATK_Variants_batch.py ${MANIFEST_PATH} /mnt $BATCH_ID

    if [ $? -ne 0 ]; then
        echo "Non-zero return code received from GATK_variants_batch.py"
        exit 1
    fi

    echo "GATK Variants pipeline completed successfully"


    rm -f /mnt/temp/*.vcf
    rm -f /mnt/temp/*.vcf.idx
    rm -f /mnt/temp/*.table
    rm -f /mnt/temp/*.table.idx
    rm -f /mnt/temp/*.tranches
    rm -f /mnt/temp/*.rscript
    rm -f /mnt/temp/*.pdf
    rm -f split_vcf_done.flag
fi
######
# Run FreeBayes script
# Similar to above, but simplified for FreeBayes running
if [ ! -n "$SKIP_FB" ]; then

    python ${SCRIPT_DIR}/FreeBayes_Variants.py ${MANIFEST_PATH} /mnt $BATCH_ID

    if [ $? -ne 0 ]; then
        echo "Non-zero return code received from FreeBayes_Variants.py"
        exit 1
    fi

    echo "FreeBayes Variants pipeline completed successfully"

    rm -f /mnt/temp/*.vcf

fi
echo "SNV/Indel variant calling completed. (batch ${BATCH_ID})"