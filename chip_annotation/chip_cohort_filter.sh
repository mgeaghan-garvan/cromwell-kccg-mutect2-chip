#!/bin/bash

PREVALENCE_THRESHOLD=${1}
CHIP_CONTAINER=${2}
# Remaining arguments are input VCFs
shift
shift
INPUT_VCFS=${@}

# Check input VCFs exist
for INPUT_VCF in ${INPUT_VCFS}; do
    if [ ! -f "${INPUT_VCF}" ]; then
        echo "Input VCF ${INPUT_VCF} does not exist"
        exit 1
    fi
done

TEMP_DIR=$(mktemp -d)
OUT_DIR=output
mkdir -p output

set -euo pipefail

# === STEP 0: Check if docker or singularity is installed ===
if [ -x "$(command -v docker)" ]; then
    echo "Docker is installed"
    CHIP_DOCKER_CMD="docker run --rm -v ${PWD}:${PWD} -v ${TEMP_DIR}:${TEMP_DIR} -w ${PWD} ${CHIP_CONTAINER}"
elif [ -x "$(command -v singularity)" ]; then
    echo "Singularity is installed"
    CHIP_DOCKER_CMD="singularity exec --containall -B ${PWD}:${PWD} -B ${TEMP_DIR}:${TEMP_DIR} --pwd ${PWD} ${CHIP_CONTAINER}"
else
    echo "Neither docker nor singularity is installed"
    exit 1
fi

# === STEP 1: Strip FILTER, INFO, and non-GT fields from input VCFs ===
for INPUT_VCF in ${INPUT_VCFS}; do
    INPUT_VCF_BN="$(basename $(basename ${INPUT_VCF} .gz) .vcf)"
    STRIPPED_VCF="${TEMP_DIR}/${INPUT_VCF_BN}.stripped.vcf.gz"
    bcftools annotate -x FILTER,INFO,FMT -O z -o ${STRIPPED_VCF} ${INPUT_VCF} 
    tabix -s 1 -b 2 -e 2 ${STRIPPED_VCF}
done

# === STEP 2: Merge input VCFs and split multi-allelic sites ===
MERGED_STRIPPED_SPLIT_VCF="${TEMP_DIR}/cohort.stripped.split.vcf"
bcftools merge -m all ${TEMP_DIR}/*.stripped.vcf.gz | \
    bcftools norm -m -any -o ${MERGED_STRIPPED_SPLIT_VCF}

# === STEP 3: Run cohort-wide filter ===
${CHIP_DOCKER_CMD} annotate_chip_cohort \
    --input_vcf ${MERGED_STRIPPED_SPLIT_VCF} \
    --prevalence_threshold ${PREVALENCE_THRESHOLD}

ANNOTATION_VCF="${TEMP_DIR}/$(basename ${MERGED_STRIPPED_SPLIT_VCF} .vcf).annotated.vcf.gz"

bcftools sort ${ANNOTATION_VCF} > ${TEMP_DIR}/cohort.annotated.sorted.vcf
bgzip -c ${TEMP_DIR}/cohort.annotated.sorted.vcf > ${TEMP_DIR}/cohort.annotated.sorted.vcf.gz
tabix -s 1 -b 2 -e 2 ${TEMP_DIR}/cohort.annotated.sorted.vcf.gz

bcftools norm -m +any -O z -o ${TEMP_DIR}/cohort.annotated.merged.sorted.vcf.gz ${TEMP_DIR}/cohort.annotated.sorted.vcf.gz
tabix -s 1 -b 2 -e 2 ${TEMP_DIR}/cohort.annotated.merged.sorted.vcf.gz

ANNOTATION_VCF="${TEMP_DIR}/cohort.annotated.merged.sorted.vcf.gz"

# === STEP 4: Annotate input VCFs ===
for INPUT_VCF in ${INPUT_VCFS}; do
    INPUT_VCF_BN="$(basename $(basename ${INPUT_VCF} .gz) .vcf)"
    bcftools annotate \
        -a ${ANNOTATION_VCF} \
        -c "=FILTER,+INFO" \
        -O z \
        -o ${OUT_DIR}/${INPUT_VCF_BN}.chip.cohort_filter.vcf.gz
    tabix -s 1 -b 2 -e 2 ${OUT_DIR}/${INPUT_VCF_BN}.chip.cohort_filter.vcf.gz
done
