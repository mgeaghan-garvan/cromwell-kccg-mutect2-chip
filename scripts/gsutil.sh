#!/bin/bash

GSUTIL_ARGS="${@}"

set -euo pipefail

if [ -n "${GSUTIL_ARGS}" ]
then
    GSUTIL_CMD="gsutil ${GSUTIL_ARGS}"
    echo "Running: ${GSUTIL_CMD}"
    eval "${GSUTIL_CMD}"
fi