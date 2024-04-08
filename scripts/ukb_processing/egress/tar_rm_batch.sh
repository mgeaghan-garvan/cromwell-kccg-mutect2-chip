#!/bin/bash
#PBS -q normal
#PBS -P np30
#PBS -l mem=8G
#PBS -l walltime=05:00:00
#PBS -l storage=gdata/tn36

# Check variables are specified
if [ -z "${BATCH}" ]; then echo "ERROR: MISSING INPUT VARIABLES"; exit 1; fi

set -euo pipefail

tar -cf ${BATCH}.tar ${BATCH}/

# Check MD5s
DIRMD5=$(tar -cf - ${BATCH}/ | md5sum | cut -d ' ' -f 1)
TARMD5=$(md5sum ${BATCH}.tar | cut -d ' ' -f 1)

if [ "${DIRMD5}" == "${TARMD5}" ]
then
  rm -r ${BATCH}/
else
  echo "BATCH ${BATCH} MD5s DONT MATCH"
  exit 1
fi

echo DONE
