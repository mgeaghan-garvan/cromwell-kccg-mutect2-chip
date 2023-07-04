#!/bin/bash

MODE=""
WORKFLOW_ID=""
CROMWELL_HOST="localhost"
CROMWELL_PORT="8007"
WDL_PATH=""
INPUT_JSON=""
BATCH_INPUT_TSV=""
OPTIONS_JSON=""
DX_WORKFLOW=""
DX_OUTPUT=""
DX_PRIORITY="low"

RUN_ID="$(date +%y%m%d_%H%M%S)"
OUTPUT_PATH="run_history/run.${RUN_ID}.txt"
mkdir -p run_history

# Usage
helpmsg() {
    echo "Manage the Cromwell server."
    echo -e "\nUsage: ./manage.sh [MODE] [--id WORKFLOW_ID] [--host CROMWELL_HOST] [--port CROMWELL_PORT] [--wdl WDL_PATH] [--input INPUT_JSON] [--batch BATCH_INPUT_TSV] [--options OPTIONS_JSON] [--dxworkflow DX_WORKFLOW] [--dxoutput DX_OUTPUT] [--dxpriority DX_PRIORITY]"
    echo -e "Display this help message: $0 -h\n"
    echo -e "\tMODE:                       Either 'submit', 'dxsubmit', 'status', or 'abort'.           REQUIRED"
    echo -e "\t--id WORKFLOW_ID:           Workflow ID for which to get status or abort.                REQUIRED for 'status' and 'abort' modes."
    echo -e "\t--host CROMWELL_HOST:       Host on which Cromwell web service is running.               Default: '${CROMWELL_HOST}'"
    echo -e "\t--port CROMWELL_PORT:       Port on which Cromwell web service is running.               Default: '${CROMWELL_PORT}'"
    echo -e "\t--wdl WDL_PATH:             Path to WDL file describing workflow to run.                 REQUIRED for 'submit' mode"
    echo -e "\t--input INPUT_JSON:         Path to JSON file describing workflow inputs.                REQUIRED for 'submit' and 'dxsubmit' mode"
    echo -e "\t--batch BATCH_INPUT_TSV:    Path to TSV file describing inputs for batch runs.           REQUIRED for 'dxsubmit' mode; OPTIONAL for 'submit' mode; if supplied, a batch job will be submitted"
    echo -e "\t--options OPTIONS_JSON:     Path to JSON file describing workflow options.               REQUIRED for 'submit' mode"
    echo -e "\t--dxworkflow DX_WORKFLOW:   DNAnexus path to workflow.                                   REQUIRED for 'dxsubmit' mode"
    echo -e "\t--dxoutput DX_OUTPUT:       DNAnexus path to place output.                               REQUIRED for 'dxsubmit' mode"
    echo -e "\t--dxpriority DX_PRIORITY:   DNAnexus run priority; either 'low', 'normal', or 'high'.    Default: '${DX_PRIORITY}'"
}

# Arguments
# Get run mode (submit, status, or abort)
FIRSTARG="${1}"
case ${FIRSTARG} in
    submit)
        MODE="submit"
        shift
        ;;
    dxsubmit)
        MODE="dxsubmit"
        shift
        ;;
    status)
        MODE="status"
        shift
        ;;
    abort)
        MODE="abort"
        shift
        ;;
    *)
        echo "Invalid mode."
        helpmsg
        exit 0
        ;;
esac

POSITIONAL=()
while [[ $# -gt 0 ]]; do
    key="${1}"
    case ${key} in
        -h|--help)
            helpmsg
            exit 0
            ;;
        --id)
            WORKFLOW_ID="${2}"
            shift
            shift
            ;;
        --host)
            CROMWELL_HOST="${2}"
            shift
            shift
            ;;
        --port)
            CROMWELL_PORT="${2}"
            shift
            shift
            ;;
        --wdl)
            WDL_PATH="${2}"
            shift
            shift
            ;;
        --input)
            INPUT_JSON="${2}"
            shift
            shift
            ;;
        --batch)
            BATCH_INPUT_TSV="${2}"
            shift
            shift
            ;;
        --options)
            OPTIONS_JSON="${2}"
            shift
            shift
            ;;
        --dep)
            DEPENDENCIES_ZIP="${2}"
            shift
            shift
            ;;
        --dxworkflow)
            DX_WORKFLOW="${2}"
            shift
            shift
            ;;
        --dxoutput)
            DX_OUTPUT="${2}"
            shift
            shift
            ;;
        --dxpriority)
            DX_PRIORITY="${2}"
            shift
            shift
            ;;
        *)
            POSITIONAL+=("$1")
            shift
            ;;
    esac
done

set -- "${POSITIONAL[@]}"

set -euo pipefail

# Check the management mode is valid
if [ ! "${MODE}" == "submit" ] && [ ! "${MODE}" == "status" ] && [ ! "${MODE}" == "abort" ] && [ ! "${MODE}" == "dxsubmit" ]
then
    echo "Invalid mode."
    exit 1
fi

# Submit mode
if [ "${MODE}" == "submit" ]
then
    for SETTING in "${CROMWELL_HOST}" "${CROMWELL_PORT}" "${WDL_PATH}" "${INPUT_JSON}" "${OPTIONS_JSON}"
    do
        if [ -z "${SETTING}" ]
        then
            echo "ERROR: Missing required argument."
            exit 1
        fi
        case "${SETTING}" in
            *" "*)
                echo "ERROR: Invalid argument: '${SETTING}'"
                exit 1
                ;;
            *","*)
                echo "ERROR: Invalid argument: '${SETTING}'"
                exit 1
                ;;
        esac
    done
    BATCH_MODE="FALSE"
    API_ENDPOINT="api/workflows/v1"
    INPUT_JSON_FILE=${INPUT_JSON}
    if [ -n "${BATCH_INPUT_TSV}" ] && [ -s "${BATCH_INPUT_TSV}" ]
    then
        BATCH_MODE="TRUE"
        API_ENDPOINT="${API_ENDPOINT}/batch"
        # Create multiple JSON files, one per sample in the batch TSV file
        BATCH_DIR="input/config/batch/${RUN_ID}"
        rm -rf ${BATCH_DIR}
        mkdir -p ${BATCH_DIR}
        python3 scripts/generate_batch_inputs.py -j ${INPUT_JSON} -b ${BATCH_INPUT_TSV} -o ${BATCH_DIR} -f singlefile
        INPUT_JSON_FILE="${BATCH_DIR}/$(basename ${INPUT_JSON} .json).batch.json"
    fi

    # Zip up dependencies
    DEPENDENCIES_ZIP="workflow/modules.zip"
    rm -f ${DEPENDENCIES_ZIP}
    cd $(dirname ${DEPENDENCIES_ZIP})
    zip $(basename ${DEPENDENCIES_ZIP}) *.wdl
    cd -
    DEPENDENCIES_STR="-F workflowDependencies=@${DEPENDENCIES_ZIP}"

    # Submit the workflow
    echo "Submitting workflow..." > ${OUTPUT_PATH}
    echo "Run ID: ${RUN_ID}" >> ${OUTPUT_PATH}
    echo "Host: ${CROMWELL_HOST}" >> ${OUTPUT_PATH}
    echo "Port: ${CROMWELL_PORT}" >> ${OUTPUT_PATH}
    echo "WDL: ${WDL_PATH}" >> ${OUTPUT_PATH}
    echo "Inputs: ${INPUT_JSON_FILE}" >> ${OUTPUT_PATH}
    echo "Workflow options: ${OPTIONS_JSON}" >> ${OUTPUT_PATH}
    CMD="curl -X POST \"http://${CROMWELL_HOST}:${CROMWELL_PORT}/${API_ENDPOINT}\"
        -H \"accept: application/json\"
        -H \"Content-Type: multipart/form-data\"
        -F \"workflowSource=@${WDL_PATH}\"
        ${DEPENDENCIES_STR}
        -F \"workflowInputs=@${INPUT_JSON_FILE}\"
        -F \"workflowOptions=@${OPTIONS_JSON}\""
    echo "Submission command: ${CMD}" >> ${OUTPUT_PATH}
    SUBRESPONSE=$(eval ${CMD})
    echo "Submission status: ${SUBRESPONSE}" >> ${OUTPUT_PATH}
    CROMWELL_ID=$(echo ${SUBRESPONSE} | grep -oP "(?<=\"id\":\")[^\"]+(?=\")")
    echo "Cromwell submission ID: ${CROMWELL_ID}" >> ${OUTPUT_PATH}
    cat ${OUTPUT_PATH}
    echo "Run details file: ${OUTPUT_PATH}"
fi

if [ "${MODE}" == "dxsubmit" ]
then
    for SETTING in "${INPUT_JSON}" "${BATCH_INPUT_TSV}" "${DX_WORKFLOW}" "${DX_OUTPUT}" "${DX_PRIORITY}"
    do
        if [ -z "${SETTING}" ]
        then
            echo "ERROR: Missing required argument."
            exit 1
        fi
        case "${SETTING}" in
            *" "*)
                echo "ERROR: Invalid argument: '${SETTING}'"
                exit 1
                ;;
            *","*)
                echo "ERROR: Invalid argument: '${SETTING}'"
                exit 1
                ;;
        esac
    done
    # Create multiple DNAnexus submission scripts, one per sample in the batch TSV file
    BATCH_DIR="input/config/dx_batch/${RUN_ID}"
    rm -rf ${BATCH_DIR}
    mkdir -p ${BATCH_DIR}
    python3 scripts/generate_batch_inputs.py -j ${INPUT_JSON} -b ${BATCH_INPUT_TSV} -o ${BATCH_DIR} -f dx --dx_workflow ${DX_WORKFLOW} --dx_destination ${DX_OUTPUT} --dx_priority ${DX_PRIORITY}
    BATCH_SCRIPTS="${BATCH_DIR}/$(basename ${INPUT_JSON} .json).*.dx.sh"

    # Submit the workflow
    echo "Submitting DNAnexus workflow..." > ${OUTPUT_PATH}
    echo "Run ID: ${RUN_ID}" >> ${OUTPUT_PATH}
    echo "Workflow: ${DX_WORKFLOW}" >> ${OUTPUT_PATH}
    echo "Submission scripts: ${BATCH_SCRIPTS}" >> ${OUTPUT_PATH}
    echo "Priority: ${DX_PRIORITY}" >> ${OUTPUT_PATH}
    echo "Output: ${DX_OUTPUT}" >> ${OUTPUT_PATH}
    echo "Submission commands:" >> ${OUTPUT_PATH}
    for BATCH_SCRIPT in ${BATCH_SCRIPTS}
    do
        CMD="bash ${BATCH_SCRIPT}"
        echo ${CMD} >> ${OUTPUT_PATH}
        eval ${CMD} &
    done
    wait
    cat ${OUTPUT_PATH}
    echo "Run details file: ${OUTPUT_PATH}"
fi

# Status or abort modes
if [ "${MODE}" == "status" ] || [ "${MODE}" == "abort" ]
then
    for SETTING in "${CROMWELL_HOST}" "${CROMWELL_PORT}" "${WORKFLOW_ID}"
    do
        if [ -z "${SETTING}" ]
        then
            echo "ERROR: Missing required argument."
            exit 1
        fi
        case "${SETTING}" in
            *" "*)
                echo "ERROR: Invalid argument: '${SETTING}'"
                exit 1
                ;;
            *","*)
                echo "ERROR: Invalid argument: '${SETTING}'"
                exit 1
                ;;
        esac
    done

    # Status mode
    if [ "${MODE}" == "status" ]
    then
        curl -X GET "http://${CROMWELL_HOST}:${CROMWELL_PORT}/api/workflows/v1/${WORKFLOW_ID}/status" \
            -H "accept: application/json"
    fi

    # Abort mode
    if [ "${MODE}" == "abort" ]
    then
        curl -X POST "http://${CROMWELL_HOST}:${CROMWELL_PORT}/api/workflows/v1/${WORKFLOW_ID}/abort" \
            -H "accept: application/json"
    fi
fi

echo ""