#!/bin/bash

# test if the first argument is unset, if not, set it to the FAMILY_ID variable
if [ ! -z "$1" ]; 
then
    FAMILY_ID=$1
    echo "FAMILY_ID is ${FAMILY_ID}"
fi

SCRIPT_DIR="01_mapping_pipeline"
MANIFEST_PATH="/path/to/manifest.txt"

PROJECT_ID="default_project"
HOSTNAME=`hostname`

# S3 bucket name for storing data
S3_BUCKET="s3://your-bucket-name"
# S3 path for data
S3_DATA_PATH="data"

LOCAL_TMP=/mnt

########
# STEP 0: Check if family already proccessed!
# - Can be overriden with "$FORCE" environmental variable (e.g., qsub -v FORCE=1)
if [ -n "$(s3cmd ls ${S3_BUCKET}/${S3_DATA_PATH}/${FAMILY_ID}/)" ]; then
    if [ ! -n "$FORCE" ]; then
        echo "Already completed processing for ${FAMILY_ID}"
        exit 0
    fi
fi

echo "Starting driver.sh for ${FAMILY_ID}"

########
# STEP 1: Init_Pipeline
# - rsync reference files
# - download input files from s3
# TODO: don't download files if checkpoints exist!
python ${SCRIPT_DIR}/1_Init_Pipeline.py ${RUFFUS_ARGS} ${MANIFEST_PATH} ${FAMILY_ID} $LOCAL_TMP $MAP_SE
CODE=$?
if [ $CODE -ne 0 ]; then
    echo "Non-zero return code ${CODE} received from 1_Init_Pipeline.py, skipping ${FAMILY_ID} and continuing on."
    rm -rf ${LOCAL_TMP}/temp/*
    rm -rf ${LOCAL_TMP}/out/*
    exit 1
fi

#######
# STEP 2: BWA_Pipeline
# - map using BWA(-mem)
# - Picard CleanSam
# - Picard FixMateInformation
# - Picard MarkDuplicates
# - Run QC Metrics
python ${SCRIPT_DIR}/2_BWA_Pipeline.py ${RUFFUS_ARGS} ${MANIFEST_PATH} ${FAMILY_ID} ${LOCAL_TMP} --bwa-mem --bwa-seed 17 $MAP_SE $PREFIX_RG
if [ $? -ne 0 ]; then
    echo "Non-zero return code received from 2_BWA_Pipeline.py"
    exit 1
fi

#######
# STEP 3: GATK Cleanup
# - Indel Realignment (family-based)
# - Base recalibration
python ${SCRIPT_DIR}/3_GATK_Pipeline.py ${RUFFUS_ARGS} ${MANIFEST_PATH} ${FAMILY_ID} ${LOCAL_TMP}
if [ $? -ne 0 ]; then
    echo "Non-zero return code received from 3_GATK_Pipeline.py"
    exit 1
fi

#######
# STEP 4: QC Pipeline
# - run qplot on all bam files and generate plots

python ${SCRIPT_DIR}/QC_Pipeline.py ${RUFFUS_ARGS} ${MANIFEST_PATH} ${FAMILY_ID} ${LOCAL_TMP}
if [ $? -ne 0 ]; then
    echo "Non-zero return code received from 4_QC_Pipeline.py"
    exit 1
fi

#######
# STEP 5: Cleanup and Upload
# - Upload to s3
# - Cleanup directories
python ${SCRIPT_DIR}/Cleanup_Pipeline.py ${RUFFUS_ARGS} ${MANIFEST_PATH} ${FAMILY_ID} ${LOCAL_TMP}
if [ $? -ne 0 ]; then
    echo "Non-zero return code received from 5_Cleanup_Pipeline.py"
    exit 1
fi

echo "Completed processing for ${FAMILY_ID}"
