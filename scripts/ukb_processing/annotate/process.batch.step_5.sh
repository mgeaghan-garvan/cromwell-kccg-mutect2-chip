#!/bin/bash
#PBS -q normal
#PBS -P tn36
#PBS -l ncpus=2
#PBS -l mem=8G
#PBS -l walltime=12:00:00
#PBS -l storage=gdata/tn36+scratch/tn36

# This script takes a batch number and run number and
# extracts the corresponding VCFs,
# creates split-multi-allelic VCFs,
# creates a sites-only VCF,
# annotates the sites-only VCF,
# imports the annotations into the split-multi-allelic VCFs, and
# subsets the VCFs for UBA1

# Check variables are specified
if [ -z "${BATCH}" ] || [ -z "${TEMPDIR}" ] || [ -z "${PROCESSDIR}" ]; then echo "ERROR: MISSING INPUT VARIABLES"; exit 1; fi

# Ensure bcftools and annovar are in the PATH variable
export PATH=${PATH}:${HOME}/software/htslib/bin/bin
export PATH=${PATH}:${HOME}/software/annovar

set -euo pipefail

cd ${PROCESSDIR}

# === STEP 5 - ARCHIVE AND REMOVE VCFs ===
SPLITVCFANNOTBATCHDIR=${TEMPDIR}/vcf/annotated/batch_${BATCH}
mkdir -p vcf
TARFILE=vcf/batch_${BATCH}.annotated_vcfs.tar
tar -cf ${TARFILE} -C ${SPLITVCFANNOTBATCHDIR} .

# Remove temporary files
rm -rf ${TEMPDIR}

echo DONE
