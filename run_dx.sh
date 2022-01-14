#!/bin/bash

# Usage
helpmsg() {
    echo "Run the Mutect2 pipeline on DNAnexus."
    echo -e "\nUsage: $0 [-w|--workflow PATH_TO_WORKFLOW] [-i|--input INPUT_JSON] [-o|--output OUTPUT_PATH]"
    echo -e "Display this help message: $0 -h\n"
    echo -e "\tPATH_TO_WORKFLOW:  Path to the workflow on the DNAnexus platform."
    echo -e "\tINPUT_JSON:        Input JSON file (use the templates in the input/ directory)."
    echo -e "\tOUTPUT_PATH:       Path to the output directory on the DNAnexus platform."
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

# Create DNAnexus input JSON file
python3 generate_dx_run_cmd.py -j ${INPUT_JSON} -d ${OUTPUT_PATH} -w ${PATH_TO_WORKFLOW} -o _dx_run.sh

echo "DNAnexus run command:"
cat _dx_run.sh
chmod +x _dx_run.sh

if [ "${DRYRUN}" == "FALSE" ]
then
    ./_dx_run.sh
fi
