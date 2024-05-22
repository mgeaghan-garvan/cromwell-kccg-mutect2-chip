#!/bin/bash

set -euo pipefail

# Assumes that analyses have been organised on DNAnexus like so:
# /DXROOT/BATCH/run_xx
# Batches are based on the first two digits of the sample IDs
# Each run is a sub-batch of 500 samples

# Check variables are specified
if [ -z "${DXROOT}" ] || [ -z "${BATCH}" ]; then echo "ERROR: MISSING INPUT VARIABLES"; exit 1; fi

dx ls "${DXROOT}/${BATCH}" | sed -e 's|/$||g'
