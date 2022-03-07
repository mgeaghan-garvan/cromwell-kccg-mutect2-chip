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
if [[ -n $1 ]]; then
        echo "Remaining arguments: "
        echo "$@"
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
sed -i -e "s/CROMWELL_PORT_TO_SED/${CROMPORT}/g" ./run.sh
sed -i -e "s/CROMWELL_PORT_TO_SED/${CROMPORT}/g" ./run_chip_only.sh
sed -i -e "s/CROMWELL_PORT_TO_SED/${CROMPORT}/g" ./run_vep_only.sh
sed -i -e "s/CROMWELL_PORT_TO_SED/${CROMPORT}/g" ./run_annovar_only.sh
sed -i -e "s/CROMWELL_PORT_TO_SED/${CROMPORT}/g" ./create_pon.sh
sed -i -e "s/CROMWELL_PORT_TO_SED/${CROMPORT}/g" ./abort.sh
if [ "${PLATFORM}" == "GCP" ]
then
    sed -i -e "s/options\.json/options\.google\.json/g" ./run.sh
    sed -i -e "s/options\.json/options\.google\.json/g" ./run_chip_only.sh
    sed -i -e "s/options\.json/options\.google\.json/g" ./run_vep_only.sh
    sed -i -e "s/options\.json/options\.google\.json/g" ./run_annovar_only.sh
    sed -i -e "s/options\.json/options\.google\.json/g" ./create_pon.sh
fi
sed -i -e "s/MULTI_TO_SED/${MULTI}/g" ./run.sh
sed -i -e "s/MULTI_TO_SED/${MULTI}/g" ./run_chip_only.sh
sed -i -e "s/MULTI_TO_SED/${MULTI}/g" ./run_vep_only.sh
sed -i -e "s/MULTI_TO_SED/${MULTI}/g" ./run_annovar_only.sh
cd workflow
zip cromwell-kccg-mutect2.multi.dep.zip cromwell-kccg-mutect2.wdl cromwell-kccg-mutect2.chip.wdl cromwell-kccg-mutect2.vep.wdl cromwell-kccg-mutect2.annovar.wdl
cd ..
