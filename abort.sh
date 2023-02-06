#!/bin/bash

RUN_ID=""
if [ -z "${1}" ] && [ -f "run_id.txt" ]
then
    RUN_ID=$(cat run_id.txt | grep -oP "(?<=\"id\":\")[^\"]+(?=\")")
else
    RUN_ID="${1}"
fi

curl -X POST "http://localhost:CROMWELL_PORT_TO_SED/api/workflows/v1/${RUN_ID}/abort"
