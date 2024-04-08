#!/bin/bash

# Check variables are specified
if [ -z "${BATCH}" ] || [ -z "${RUN}" ]; then echo "ERROR: MISSING INPUT VARIABLES"; exit 1; fi

set -euo pipefail

cd ${BATCH}

TARFILE=${BATCH}.${RUN}.tar.gz
rm ${TARFILE}

cd -

echo DONE
