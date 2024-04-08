#!/bin/bash
#PBS -q normal
#PBS -P tn36
#PBS -l ncpus=4
#PBS -l mem=8G
#PBS -l walltime=00:30:00
#PBS -l storage=gdata/tn36+scratch/tn36

# This script takes a batch number and run number and
# extracts the corresponding VCFs,
# creates split-multi-allelic VCFs,
# creates a sites-only VCF,
# annotates the sites-only VCF,
# imports the annotations into the split-multi-allelic VCFs, and
# subsets the VCFs for UBA1

# Check variables are specified
if [ -z "${BATCH}" ] || [ -z "${RUN}" ] || [ -z "${TEMPDIR}" ] || [ -z "${PROCESSDIR}" ]; then echo "ERROR: MISSING INPUT VARIABLES"; exit 1; fi

# Ensure bcftools and annovar are in the PATH variable
export PATH=${PATH}:${HOME}/software/htslib/bin/bin
export PATH=${PATH}:${HOME}/software/annovar

set -euo pipefail

cd ${PROCESSDIR}

VCFDIR=${TEMPDIR}/vcf/m2

VCFRUNDIR=${VCFDIR}/batch_${BATCH}/run_${RUN}
# === STEP 2 - CREATE SPLIT-MULTI-ALLELIC VCFs ===
SPLITVCFDIR=${TEMPDIR}/vcf/split/batch_${BATCH}/run_${RUN}
rm -rf ${SPLITVCFDIR}
mkdir -p ${SPLITVCFDIR}

SOVCFDIR=${TEMPDIR}/vcf/sites_only
mkdir -p ${SOVCFDIR}

SOVCFTEMPDIR=${SOVCFDIR}/batch_${BATCH}/run_${RUN}
rm -rf ${SOVCFTEMPDIR}
mkdir -p ${SOVCFTEMPDIR}

for VCF in ${VCFRUNDIR}/*/*.vcf.gz
do
    SAMPLE=$(basename ${VCF} .vcf.gz)
    SPLITVCF=${SPLITVCFDIR}/${SAMPLE}.split.vcf.gz
    bcftools norm --threads 4 -m -any -O z -o ${SPLITVCF} ${VCF}
    tabix -s1 -b2 -e2 ${SPLITVCF}

    # === STEP 2.1 - CREATE SITES-ONLY VCF ===
    SOVCF=${SOVCFTEMPDIR}/${SAMPLE}.so.vcf.gz
    bcftools annotate --threads 4 -x FMT,INFO,FILTER ${SPLITVCF} | cut -f 1-8 | bgzip -c > ${SOVCF}
    tabix -s1 -b2 -e2 ${SOVCF}
done

SOVCF=${SOVCFDIR}/batch_${BATCH}_run_${RUN}.so.vcf.gz
rm -f ${SOVCF}

bcftools concat --threads 4 -a -d exact ${SOVCFTEMPDIR}/*.so.vcf.gz | bgzip -c > ${SOVCF}
tabix -s1 -b2 -e2 ${SOVCF}

rm -rf ${SOVCFTEMPDIR}

echo "STEP 2 - SPLIT VCFS AND CREATE SITES ONLY VCF - DONE"
