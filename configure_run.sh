#!/bin/bash
# set -euo pipefail
# set -x

# Usage
helpmsg() {
    echo "Configure Cromwell to run the Mutect2 pipeline."
    echo -e "\nUsage: $0 [-n|--name RUN_NAME] [-p|--cromport CROMWELL_PORT] [-f|--platform PLATFORM] [-d|--dryrun] [-m|--multi]"
    echo -e "Display this help message: $0 -h\n"
    echo -e "\tRUN_NAME:      Name of the run.                                            (Default: 'run')."
    echo -e "\tCROMWELL_PORT: Port where Cromwell should run.                             (Default: '8007')"
    echo -e "\tPLATFORM:      Platform on which workflow should be run.                   (Options: 'HPC', 'GCP'. Default: 'HPC')"
    echo -e "\t[-m|--multi]:  Run in multi-sample batch mode (requires inputFiles.tsv)."
    echo -e "\t[-d|--dryrun]: Print settings to screen without making changes."
}

# Display help message if there are no arguments
# Commented out to allow for default settings
# if [ $# -le 1 ];
# then
#     helpmsg
#     exit 1
# fi

# Set defaults
# Host, port, database and platform defaults
RUNNAME="run"
CROMPORT="8007"
PLATFORM="HPC"
MULTI="FALSE"
DRYRUN="0"
DX_PATH_TO_WORKFLOW=""
DX_INPUT_JSON=""
DX_BATCH_INPUT_FILE=""
DX_OUTPUT_PATH=""

# Arguments
POSITIONAL=()
while [[ $# -gt 0 ]]; do
        key="$1"
        case ${key} in
                -h|--help)
                        helpmsg
                        exit 0
                        ;;
                -n|--name)
                        RUNNAME="$2"
                        shift
                        shift
                        ;;
                -p|--cromport)
                        CROMPORT="$2"
                        shift
                        shift
                        ;;
                -f|--platform)
                        PLATFORM="$2"
                        # Set allowed platforms
                        case ${PLATFORM} in
                            "HPC")
                                ;;
                            "GCP")
                                ;;
                            "DX")
                                ;;
                            *)
                                echo -e "ERROR: Invalid platform selection.\n"
                                helpmsg
                                exit 1
                                ;;
                        esac
                        shift
                        shift
                        ;;
                -m|--multi)
                        MULTI="TRUE"
                        shift
                        ;;
                -w|--workflow)
                        DX_PATH_TO_WORKFLOW="$2"
                        shift
                        shift
                        ;;
                -i|--input)
                        DX_INPUT_JSON="$2"
                        shift
                        shift
                        ;;
                -b|--batch)
                        DX_BATCH_INPUT_FILE="$2"
                        shift
                        shift
                        ;;
                -o|--output)
                        DX_OUTPUT_PATH="$2"
                        shift
                        shift
                        ;;
                -d|--dryrun)
                        DRYRUN="1"
                        shift
                        ;;
                *)
                        POSITIONAL+=("$1")
                        shift
                        ;;
        esac
done

set -- "${POSITIONAL[@]}"

echo "Run name            = ${RUNNAME}"
echo "Cromwell port       = ${CROMPORT}"
echo "Platform            = ${PLATFORM}"
echo "Multi-sample mode   = ${MULTI}"
if [ "${PLATFORM}" == "DX" ]; then
    echo "DNAnexus workflow   = ${DX_PATH_TO_WORKFLOW}"
    echo "DNAnexus input JSON = ${DX_INPUT_JSON}"
    echo "DNAnexus batch file = ${DX_BATCH_INPUT_FILE}"
    echo "DNAnexus output     = ${DX_OUTPUT_PATH}"
fi
if [[ -n $1 ]]; then
        echo "Remaining arguments: "
        echo "$@"
fi

if [ "${PLATFORM}" == "DX" ]
then
    # Run checks to make sure all the information is provided to run on DNAnexus
    if [ "${DX_PATH_TO_WORKFLOW}" == "" ] || [ "${DX_INPUT_JSON}" == "" ] || [ "${DX_OUTPUT_PATH}" == "" ]
    then
        echo "Missing DNAnexus configuration information. Exiting."
        exit 1
    fi
    if [ ! -f "${DX_INPUT_JSON}" ]; then echo "Invalid path to input JSON file!"; exit 1; fi
    if [ "${DX_BATCH_INPUT_FILE}" != "" ] && [ ! -f "${DX_BATCH_INPUT_FILE}" ]; then echo "Invalid path to batch TSV file!"; exit 1; fi
    # Construct command to generate DNAnexus run script(s)
    CMD=""
    mkdir -p ./scripts/dx
    if [ "${DX_BATCH_INPUT_FILE}" == "" ]
    then
        CMD="python3 ./scripts/generate_dx_run_cmd.py -j ${DX_INPUT_JSON} -d ${DX_OUTPUT_PATH} -w ${DX_PATH_TO_WORKFLOW} -o ./scripts/dx/_dx_run.sh"
    else
        CMD="python3 ./scripts/generate_dx_run_cmd.py -j ${INPUT_JSON} -b ${BATCH_INPUT_FILE} -d ${OUTPUT_PATH} -w ${PATH_TO_WORKFLOW} -o ./scripts/dx/_dx_run.sh"
    fi
    if [ "${DRYRUN}" == "1" ]
    then
            # Echo python command to generate run script(s)
            "Command to generate DNAnexus run script(s):"
            echo ${CMD}
    else
        # Create DNAnexus input JSON file(s)
        CMD_FILE_LIST=$(eval ${CMD})
        for FILE in ${CMD_FILE_LIST}
        do
            echo bash ${FILE}
        done > ./scripts/run_dx.sh
    fi
    exit 0  # The rest of the configuration script only applies to HPC and GCP runs
fi

if [ "${DRYRUN}" == "1" ]
then
        exit 0
fi

# Configure options files
# '#' delimiters are used to enable path substitution
sed -i -e "s#PWD_TO_SED#${PWD}#g" ./workflow/options.json
sed -i -e "s#PWD_TO_SED#${PWD}#g" -e "s/ID_TO_SED/${RUNNAME}/g" ./workflow/options.google.json

# Set up output directories
mkdir -p workflow_out
mkdir -p workflow_logs
mkdir -p workflow_call_logs
if [ "${PLATFORM}" == "GCP" ]
then
    mkdir -p gcp_logs
fi

# Configure run.sh
sed -i -e "s/CROMWELL_PORT_TO_SED/${CROMPORT}/g" ./scripts/run_full.sh
sed -i -e "s/CROMWELL_PORT_TO_SED/${CROMPORT}/g" ./scripts/run_chip_only.sh
sed -i -e "s/CROMWELL_PORT_TO_SED/${CROMPORT}/g" ./scripts/run_vep_only.sh
sed -i -e "s/CROMWELL_PORT_TO_SED/${CROMPORT}/g" ./scripts/run_annovar_only.sh
sed -i -e "s/CROMWELL_PORT_TO_SED/${CROMPORT}/g" ./scripts/create_pon.sh
sed -i -e "s/CROMWELL_PORT_TO_SED/${CROMPORT}/g" ./abort.sh
if [ "${PLATFORM}" == "GCP" ]
then
    sed -i -e "s/options\.json/options\.google\.json/g" ./scripts/run_full.sh
    sed -i -e "s/options\.json/options\.google\.json/g" ./scripts/run_chip_only.sh
    sed -i -e "s/options\.json/options\.google\.json/g" ./scripts/run_vep_only.sh
    sed -i -e "s/options\.json/options\.google\.json/g" ./scripts/run_annovar_only.sh
    sed -i -e "s/options\.json/options\.google\.json/g" ./scripts/create_pon.sh
fi
sed -i -e "s/MULTI_TO_SED/${MULTI}/g" ./scripts/run_full.sh
sed -i -e "s/MULTI_TO_SED/${MULTI}/g" ./scripts/run_chip_only.sh
sed -i -e "s/MULTI_TO_SED/${MULTI}/g" ./scripts/run_vep_only.sh
sed -i -e "s/MULTI_TO_SED/${MULTI}/g" ./scripts/run_annovar_only.sh
cd workflow
zip cromwell-kccg-mutect2.multi.dep.zip cromwell-kccg-mutect2.wdl cromwell-kccg-mutect2.chip.wdl cromwell-kccg-mutect2.vep.wdl cromwell-kccg-mutect2.annovar.wdl
cd ..
