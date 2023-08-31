#!/bin/bash

WORKDIR=${PWD}
DBDIR=${1}

set -euo pipefail

cd ${DBDIR}

tar -czvf ${WORKDIR}/annovar_files.tar.gz *

echo "Annovar database files archived to ${WORKDIR}/annovar_files.tar.gz"
