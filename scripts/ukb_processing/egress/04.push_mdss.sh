#!/bin/bash
#PBS -q copyq
#PBS -P np30
#PBS -l mem=8G
#PBS -l walltime=10:00:00
#PBS -l storage=gdata/tn36+massdata/tn36

# Check variables are specified
if [ -z "${BATCH}" ] || [ -z "${RUN}" ] || [ -z "${PROJ}" ] || [ -z "${DEST}" ]; then echo "ERROR: MISSING INPUT VARIABLES"; exit 1; fi

set -euo pipefail

cd ${BATCH}

TARFILE=${BATCH}.${RUN}.tar.gz

# Push to MDSS
mdss -P ${PROJ} mkdir ${DEST}
CMD="mdss -P ${PROJ} put ${TARFILE} ${DEST}/${TARFILE}"
echo ${CMD}
time eval ${CMD}

# Verify
CMD="mdss -P ${PROJ} verify ${DEST}/${TARFILE}"
echo ${CMD}
time eval ${CMD}

cd -

echo DONE

