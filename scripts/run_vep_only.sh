#!/bin/bash

MULTIMODE="MULTI_TO_SED"

if [ "${MULTIMODE}" == "TRUE" ]
then
    curl -X POST "http://CROMWELL_HOST_TO_SED:CROMWELL_PORT_TO_SED/api/workflows/v1" \
        -H "accept: application/json" \
        -H "Content-Type: multipart/form-data" \
        -F "workflowSource=@./workflow/cromwell-kccg-mutect2.vep.multi.wdl" \
        -F "workflowDependencies=@./workflow/cromwell-kccg-mutect2.multi.dep.zip" \
        -F "workflowInputs=@./input/inputs.vep.multi.json;type=application/json" \
        -F "workflowOptions=@./workflow/options.json;type=application/json" | tee run_id.txt
else
    curl -X POST "http://CROMWELL_HOST_TO_SED:CROMWELL_PORT_TO_SED/api/workflows/v1" \
        -H "accept: application/json" \
        -H "Content-Type: multipart/form-data" \
        -F "workflowSource=@./workflow/cromwell-kccg-mutect2.vep.wdl" \
        -F "workflowDependencies=@./workflow/cromwell-kccg-mutect2.multi.dep.zip" \
        -F "workflowInputs=@./input/inputs.vep.json;type=application/json" \
        -F "workflowOptions=@./workflow/options.json;type=application/json" | tee run_id.txt
fi


echo ""
echo "WORKFLOW SUBMITTED"
