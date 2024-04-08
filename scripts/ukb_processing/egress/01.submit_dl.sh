#!/bin/bash

# Check variables are specified
if [ -z "${DXROOT}" ] || [ -z "${BATCH}" ] || [ -z "${RUN}" ]; then echo "ERROR: MISSING INPUT VARIABLES"; exit 1; fi

set -euo pipefail

rm -rf ${BATCH}/${RUN} && mkdir -p ${BATCH}/${RUN}
rm -rf split/${BATCH}/${RUN} && mkdir -p split/${BATCH}/${RUN}

SAMPLESTMPFILE=".${BATCH}.${RUN}.samples.txt"
dx ls "${DXROOT}/${BATCH}/${RUN}" | sed -e 's|/$||g' > ${SAMPLESTMPFILE}
split -l 50 --numeric-suffixes=1 --suffix-length=2 --additional-suffix=.txt ${SAMPLESTMPFILE} split/${BATCH}/${RUN}/samples.split_
rm ${SAMPLESTMPFILE}

mkdir -p logs/dl/${BATCH}/${RUN}

# Loop through split files and submit jobs
for f in split/${BATCH}/${RUN}/samples.split_*.txt
do
    qsub -o logs/dl/${BATCH}/${RUN} -e logs/dl/${BATCH}/${RUN} -v SAMPLESFILE=${PWD}/$f,BATCH=${BATCH},RUN=${RUN} dl.sh
done

echo "JOBS SUBMITTED"
