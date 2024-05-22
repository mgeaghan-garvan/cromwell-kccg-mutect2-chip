#!/bin/bash

set -euo pipefail

# Check variables are specified
if [ -z "${DXROOT}" ] || [ -z "${BATCH}" ] || [ -z "${RUN}" ]; then echo "ERROR: MISSING INPUT VARIABLES"; exit 1; fi

rm -rf split_rm_dx/${BATCH}/${RUN} && mkdir -p split_rm_dx/${BATCH}/${RUN}

SAMPLESTMPFILE=".${BATCH}.${RUN}.samples.txt"
dx ls "${DXROOT}/${BATCH}/${RUN}" | sed -e 's|/$||g' > ${SAMPLESTMPFILE}
split -l 50 --numeric-suffixes=1 --suffix-length=2 --additional-suffix=.txt ${SAMPLESTMPFILE} split_rm_dx/${BATCH}/${RUN}/samples.split_
rm ${SAMPLESTMPFILE}

mkdir -p logs/rm_dx/${BATCH}/${RUN}

# Loop through split files and submit jobs
for f in split_rm_dx/${BATCH}/${RUN}/samples.split_*.txt
do
    qsub -o logs/rm_dx/${BATCH}/${RUN} -e logs/rm_dx/${BATCH}/${RUN} -v SAMPLESFILE=${PWD}/$f,BATCH=${BATCH},RUN=${RUN} rm_dx.sh
done

echo "JOBS SUBMITTED"
