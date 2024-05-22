#!/bin/bash
#PBS -q normal
#PBS -P tn36
#PBS -l ncpus=1
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

SOVCFDIR=${TEMPDIR}/vcf/sites_only
mkdir -p ${SOVCFDIR}

SOVCF=${SOVCFDIR}/batch_${BATCH}.so.vcf.gz
rm -f ${SOVCF}

bcftools concat --threads 1 -a -d exact ${SOVCFDIR}/batch_*_run_*.so.vcf.gz | bgzip -c > ${SOVCF}
tabix -s1 -b2 -e2 ${SOVCF}

# === STEP 3 - ANNOTATE SITES-ONLY VCF ===
ANNOVARDIR=${HOME}/software/annovar
ANNOVAR=${ANNOVARDIR}/table_annovar.pl
DBDIR=/g/data/tn36/ukb/annovar_db

ANNOVARTEMPDIR=${TEMPDIR}/annovar_${BATCH}
rm -rf ${ANNOVARTEMPDIR}
mkdir ${ANNOVARTEMPDIR}

PREFIX=${ANNOVARTEMPDIR}/batch_${BATCH}.annotate

perl ${ANNOVAR} \
    ${SOVCF} \
    ${DBDIR} \
    -buildver hg38 \
    -out ${PREFIX} \
    -protocol cosmic95_coding,cosmic95_noncoding,dbnsfp42a,clinvar_20221231,gnomad211_exome,refGene \
    -operation f,f,f,f,f,g \
    -nastring . \
    -vcfinput \
    -polish \
    -thread 1 \
    -maxgenethread 1

ANNOVARVCF=${ANNOVARTEMPDIR}/batch_${BATCH}.annovar.vcf.gz
rm -f ${ANNOVARVCF}

bgzip -c ${PREFIX}.hg38_multianno.vcf > ${ANNOVARVCF}
tabix -s1 -b2 -e2 ${ANNOVARVCF}

rm -rf ${SOVCFDIR}/batch_*_run_*.so.vcf.gz

echo "STEP 3 - ANNOVAR - DONE"
