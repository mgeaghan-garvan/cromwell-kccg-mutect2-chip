#!/bin/bash
#PBS -q normal
#PBS -P tn36
#PBS -l ncpus=16
#PBS -l mem=64G
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
if [ -z "${BATCH}" ] || [ -z "${RUN}" ] || [ -z "${TEMPDIR}" ] || [ -z "${PROCESSDIR}" ]; then echo "ERROR: MISSING INPUT VARIABLES"; exit 1; fi

# Ensure bcftools and annovar are in the PATH variable
export PATH=${PATH}:${HOME}/software/htslib/bin/bin
export PATH=${PATH}:${HOME}/software/annovar

set -euo pipefail

cd ${PROCESSDIR}

# === STEP 4 - IMPORT ANNOTATIONS INTO SPLIT-MULTI-ALLELIC VCFs ===
SPLITVCFBATCHDIR=${TEMPDIR}/vcf/split/batch_${BATCH}
SPLITVCFANNOTBATCHDIR=${TEMPDIR}/vcf/annotated/batch_${BATCH}
SPLITVCFDIR=${SPLITVCFBATCHDIR}/run_${RUN}
SPLITVCFANNOTDIR=${SPLITVCFANNOTBATCHDIR}/run_${RUN}
rm -rf ${SPLITVCFANNOTDIR}
mkdir -p ${SPLITVCFANNOTDIR}

ls -1 ${SPLITVCFDIR}/*.split.vcf.gz > ${TEMPDIR}/batch_${BATCH}_run_${RUN}.vcf_list.txt
cd ${TEMPDIR}
split -l 100 batch_${BATCH}_run_${RUN}.vcf_list.txt batch_${BATCH}_run_${RUN}.vcf_list.split
cd -

for VCF_LIST in ${TEMPDIR}/batch_${BATCH}_run_${RUN}.vcf_list.split*
do
    while read VCF
    do
        FINALVCF="${SPLITVCFANNOTDIR}/$(basename ${VCF} .vcf.gz).annotated.vcf.gz"
        bcftools annotate --threads 2 -a ${ANNOVARVCF} -c CHROM,POS,REF,ALT,INFO ${VCF} | bgzip -c > ${FINALVCF} &
    done < ${VCF_LIST}
    wait
done
for FINALVCF in ${SPLITVCFANNOTDIR}/*.annotated.vcf.gz
do
    tabix -s1 -b2 -e2 ${FINALVCF}
done

# === STEP 4.1 - SUBSET FOR UBA1 REGION ===
for VCF in ${SPLITVCFANNOTDIR}/*.annotated.vcf.gz
do
    UBA1VCF="${SPLITVCFANNOTDIR}/$(basename ${VCF} .vcf.gz).uba1.vcf.gz"
    bcftools view --threads 16 -r chrX:47190861-47215128 ${VCF} | bgzip -c > ${UBA1VCF}
    tabix -s1 -b2 -e2 ${UBA1VCF}
done

echo DONE
