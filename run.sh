#!/bin/bash

curl -X POST "http://localhost:CROMWELL_PORT_TO_SED/api/workflows/v1" \
    -H "accept: application/json" \
    -H "Content-Type: multipart/form-data" \
    -F "workflowSource=@./workflow/cromwell-kccg-mutect2.wdl" \
    -F "workflowInputs=@./workflow/inputs.json;type=application/json" \
    -F "workflowOptions=@./workflow/options.json;type=application/json"

echo ""
echo "WORKFLOW SUBMITTED"
