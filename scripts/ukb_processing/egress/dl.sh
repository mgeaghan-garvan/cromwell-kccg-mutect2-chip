#!/bin/bash
#PBS -q copyq
#PBS -P np30
#PBS -l mem=8G
#PBS -l walltime=01:00:00
#PBS -l storage=gdata/tn36

eval "$(~/software/miniforge-pypy3/bin/conda shell.bash hook)" 

conda activate dx

# Check variables are specified
if [ -z "${SAMPLESFILE}" ] || [ -z "${DXROOT}" ] || [ -z "${BATCH}" ] || [ -z "${RUN}" ]; then echo "ERROR: MISSING INPUT VARIABLES"; exit 1; fi

set -euo pipefail

cd ${BATCH}/${RUN}

while read SAMPLE
do
    dx download -r --lightweight "${DXROOT}/${BATCH}/${RUN}/${SAMPLE}"
done < ${SAMPLESFILE}

cd -

echo DONE
