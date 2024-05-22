#!/bin/bash
#PBS -q normal
#PBS -P np30
#PBS -l mem=8G
#PBS -l walltime=01:00:00
#PBS -l storage=gdata/tn36

# Check variables are specified
if [ -z "${BATCH}" ] || [ -z "${RUN}" ]; then echo "ERROR: MISSING INPUT VARIABLES"; exit 1; fi

set -euo pipefail

cd ${BATCH}/${RUN}

tar -czf ../${BATCH}.${RUN}.tar.gz .

cd -

echo DONE
