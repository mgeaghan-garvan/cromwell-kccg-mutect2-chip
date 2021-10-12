#!/bin/bash

curl -X POST "http://localhost:CROMWELL_PORT_TO_SED/api/workflows/v1" \
    -H "accept: application/json" \
    -H "Content-Type: multipart/form-data" \
    -F "workflowSource=@./workflow/cromwell-kccg-mutect2.pon.wdl" \
    -F "workflowDependencies=@./workflow/cromwell-kccg-mutect2.multi.dep.zip" \
    -F "workflowInputs=@./workflow/inputs.pon.json;type=application/json" \
    -F "workflowOptions=@./workflow/options.json;type=application/json" | tee run_id.txt

echo ""
echo "WORKFLOW SUBMITTED"
