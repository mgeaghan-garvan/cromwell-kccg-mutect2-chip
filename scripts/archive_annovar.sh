set -euo pipefail

# Usage
helpmsg() {
    echo "Archive annovar scripts and reference files for use in the CHIP pipeline."
    echo -e "\nUsage: $0 [-s|--scripts SCRIPT_DIR] [-a|--annotations ANNOTATION_DIR]"
    echo -e "Display this help message: $0 -h\n"
    echo -e "\tSCRIPT_DIR:      Directory path containing the annovar perl scripts."
    echo -e "\tANNOTATION_DIR:  Directory path containing the annovar annotation files."
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
SCRIPT_DIR="run"
ANNOTATION_DIR="8007"

# Arguments
POSITIONAL=()
while [[ $# -gt 0 ]]; do
        key="$1"
        case ${key} in
                -h|--help)
                        helpmsg
                        exit 0
                        ;;
                -s|--scripts)
                        SCRIPT_DIR="$2"
                        shift
                        shift
                        ;;
                -a|--annotations)
                        ANNOTATION_DIR="$2"
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

mkdir annovar_files

cp ${SCRIPT_DIR}/*.pl annovar_files/
cp -r ${ANNOTATION_DIR}/* annovar_files

tar -czvf annovar_files.tar.gz annovar_files/
rm -rf annovar_files