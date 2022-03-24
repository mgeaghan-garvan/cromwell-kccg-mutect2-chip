#!/bin/bash

# Usage
helpmsg() {
    echo "Run the Mutect2 pipeline on DNAnexus."
    echo -e "\nUsage: $0 [-w|--workflow PATH_TO_WORKFLOW] [-i|--input INPUT_JSON] [-b|--batch BATCH_INPUT_FILE] [-o|--output OUTPUT_PATH]"
    echo -e "Display this help message: $0 -h\n"
    echo -e "\tPATH_TO_WORKFLOW:  Path to the workflow on the DNAnexus platform."
    echo -e "\tINPUT_JSON:        Input JSON file (use the templates in the input/ directory)."
    echo -e "\tOUTPUT_PATH:       Path to the output directory on the DNAnexus platform."
    echo -e "\tBATCH_INPUT_FILE:  Path to a batch input TSV file for running multiple single-sample workflows in parallel. See the README and the template file at ./input/inputFiles.tsv for information on the required format."
    echo -e "\t[-d|--dryrun]:     Perform a dry-run by printing the DNAnexus run command to the terminal without running."
}

# Display help message if there are no arguments
if [ $# -le 1 ];
then
    helpmsg
    exit 1
fi

# Set defaults
DRYRUN="FALSE"
PATH_TO_WORKFLOW=""
INPUT_JSON=""
BATCH_INPUT_FILE=""
OUTPUT_PATH=""

# Arguments
POSITIONAL=()
while [[ $# -gt 0 ]]; do
        key="$1"
        case ${key} in
                -h|--help)
                        helpmsg
                        exit 0
                        ;;
                -w|--workflow)
                        PATH_TO_WORKFLOW="$2"
                        shift
                        shift
                        ;;
                -i|--input)
                        INPUT_JSON="$2"
                        shift
                        shift
                        ;;
                -b|--batch)
                        BATCH_INPUT_FILE="$2"
                        shift
                        shift
                        ;;
                -o|--output)
                        OUTPUT_PATH="$2"
                        shift
                        shift
                        ;;
                -d|--dryrun)
                        DRYRUN="TRUE"
                        shift
                        ;;
                *)
                        POSITIONAL+=("$1")
                        shift
                        ;;
        esac
done

set -- "${POSITIONAL[@]}"

if [ "${PATH_TO_WORKFLOW}" == "" ] || [ "${INPUT_JSON}" == "" ] || [ "${OUTPUT_PATH}" == "" ]; then echo "Missing parameters!"; exit; fi
if [ ! -f "${INPUT_JSON}" ]; then echo "Invalid path to input JSON file!"; exit; fi
if [ "${BATCH_INPUT_FILE}" != "" ] && [ ! -f "${BATCH_INPUT_FILE}" ]; then echo "Invalid path to batch TSV file!"; exit; fi

# Create DNAnexus input JSON file
if [ "${BATCH_INPUT_FILE}" == "" ]
then
    python3 generate_dx_run_cmd.py -j ${INPUT_JSON} -d ${OUTPUT_PATH} -w ${PATH_TO_WORKFLOW} -o _dx_run.sh
    echo "DNAnexus run command:"
    cat _dx_run.sh
    chmod +x _dx_run.sh
    if [ "${DRYRUN}" == "FALSE" ]
    then
        ./_dx_run.sh
    fi
else
    CMD_FILE_LIST=$(python3 generate_dx_run_cmd.py -j ${INPUT_JSON} -b ${BATCH_INPUT_FILE} -d ${OUTPUT_PATH} -w ${PATH_TO_WORKFLOW} -o _dx_run.sh)
    for FILE in ${CMD_FILE_LIST}
    do
        echo "===== ${FILE} ====="
        echo "DNAnexus run command:"
        cat ${FILE}
        if [ "${DRYRUN}" == "FALSE" ]
        then
            bash ${FILE}
        fi
        echo "=========="
    done
fi

