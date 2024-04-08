#!/bin/bash
#PBS -q normal
#PBS -P np30
#PBS -l mem=8G
#PBS -l walltime=01:00:00
#PBS -l storage=gdata/tn36

# Check variables are specified
if [ -z "${BATCH}" ] || [ -z "${RUN}" ]; then echo "ERROR: MISSING INPUT VARIABLES"; exit 1; fi

set -euo pipefail

cd ${BATCH}

rm -rf ${RUN}_tmp && mv ${RUN} ${RUN}_tmp
mkdir -p ${RUN}

for VCF in ${RUN}_tmp/*/*-filtered.vcf.gz*
do
    SAMPLE="$(basename $(dirname ${VCF}))"
    mkdir -p ${RUN}/${SAMPLE}
    mv ${VCF} ${RUN}/${SAMPLE}/
done

rm -r ${RUN}_tmp

cd -

echo DONE
