#!/bin/bash

# Usage
helpmsg() {
    echo "Configure Cromwell to run the Mutect2 pipeline."
    echo -e "\nUsage: $0 [-m|--mode RUN_MODE]"
    echo -e "Display this help message: $0 -h\n"
    echo -e "\tRUN_MODE:      Run mode - one of: 'pon', 'full', 'chip', 'annovar', 'vep', 'dx'. (Default: 'full')."
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
RUNMODE="full"

# Arguments
POSITIONAL=()
while [[ $# -gt 0 ]]; do
        key="$1"
        case ${key} in
                -h|--help)
                        helpmsg
                        exit 0
                        ;;
                -m|--mode)
                        RUNMODE="$2"
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

echo "Run mode            = ${RUNMODE}"

if [[ -n $1 ]]; then
        echo "Remaining arguments: "
        echo "$@"
fi

if [ "${RUNMODE}" == "" ]
then
    echo "Invalid run mode selected!"
    exit 1
fi

if [ "${RUNMODE}" == "full" ]
then
    bash ./scripts/run_full.sh
    exit 0
fi

if [ "${RUNMODE}" == "pon" ]
then
    bash ./scripts/create_pon.sh
    exit 0
fi

if [ "${RUNMODE}" == "chip" ]
then
    bash ./scripts/run_chip_only.sh
    exit 0
fi

if [ "${RUNMODE}" == "annovar" ]
then
    bash ./scripts/run_annovar_only.sh
    exit 0
fi

if [ "${RUNMODE}" == "vep" ]
then
    bash ./scripts/run_vep_only.sh
    exit 0
fi

if [ "${RUNMODE}" == "dx" ]
then
    bash ./scripts/run_dx.sh
    exit 0
fi

# If run mode was valid, the script should exit before now
echo "Invalid run mode selected!"
exit 1
