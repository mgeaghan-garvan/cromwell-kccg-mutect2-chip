#!/bin/bash
#PBS -q normal
#PBS -P tn36
#PBS -l ncpus=2
#PBS -l mem=8G
#PBS -l walltime=00:10:00
#PBS -l storage=gdata/tn36+scratch/tn36

# This script takes a batch number and run number and
# extracts the corresponding VCFs,
# creates split-multi-allelic VCFs,
# creates a sites-only VCF,
# annotates the sites-only VCF,
# imports the annotations into the split-multi-allelic VCFs, and
# subsets the VCFs for UBA1

# Check variables are specified
if [ -z "${INPUTDIR}" ] || [ -z "${PROCESSDIR}" ] || [ -z "${BATCH}" ] || [ -z "${PROJ}" ]; then echo "ERROR: MISSING INPUT VARIABLES"; exit 1; fi

# Ensure bcftools and annovar are in the PATH variable
export PATH=${PATH}:${HOME}/software/htslib/bin/bin
export PATH=${PATH}:${HOME}/software/annovar

# Make a temporary directory in /scratch/tn36 and set TMPDIR
TMPDIR=$(mktemp -d --tmpdir=/scratch/${PROJ}/)
echo ${TMPDIR}

set -euo pipefail

cd ${PROCESSDIR}

# === STEP 1 - EXTRACT VCFs ===

TARFILE=${INPUTDIR}/batch_${BATCH}.tar

VCFDIR=${TMPDIR}/vcf/m2
rm -rf ${VCFDIR}
mkdir -p ${VCFDIR}

tar -xf ${TARFILE} -C ${VCFDIR} batch_${BATCH}/

# Echo all runs in batch
for RUNDIR in ${VCFDIR}/batch_${BATCH}/run_*
do
    basename ${RUNDIR}
done

echo "STEP 1 - EXTRACTION - DONE"
