#!/bin/bash
# set -euo pipefail
# set -x

# Usage
helpmsg() {
    echo "Configure Cromwell to run the Mutect2 pipeline."
    echo -e "\nUsage: $0 [-H|--dbhost DB_HOSTNAME] [-P|--dbport DB_PORT] [-n|--dbname DB_NAME] [-p|--cromport CROMWELL_PORT] [-f|--platform PLATFORM] [-M|--mysql MYSQL_DIR] [-R|--mysql_rundir MYSQL_RUN_DIR] [-c|--cromwell CROMWELL_JAR] [-d|--dryrun]"
    echo -e "Display this help message: $0 -h\n"
    echo -e "\tDB_HOSTNAME:   Hostname where MySQL server is running.            (Default: '0.0.0.0')"
    echo -e "\tDB_PORT:       Port for MySQL server on host.                     (Default: '40008')"
    echo -e "\tDB_NAME:       Name for MySQL database.                           (Default: 'cromwell')."
    echo -e "\tCROMWELL_PORT: Port where Cromwell should run.                    (Default: '8007')"
    echo -e "\tPLATFORM:      Platform on which workflow should be run.          (Options: 'HPC', 'GCP'. Default: 'HPC')"
    echo -e "\tMYSQL_DIR:     Root directory for MySQL installation.             (Default: '/home/glsai/mysql/mysql-5.7.27-linux-glibc2.12-x86_64/')"
    echo -e "\tMYSQL_RUN_DIR: Run directory for MySQL.                           (Default: '/home/glsai/mysql/mysql-5.7.27-linux-glibc2.12-x86_64/')"
    echo -e "\tCROMWELL_JAR:  Location of Cromwell JAR file.                     (Default: '/share/ClusterShare/software/contrib/micgea/cromwell/68.1/cromwell-68.1.jar')"
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
# MySQL location
MYSQL=/home/glsai/mysql/mysql-5.7.27-linux-glibc2.12-x86_64/
MYSQL_RUNDIR=/home/glsai/mysql/mysql-5.7.27-linux-glibc2.12-x86_64/
# Cromwell location
# CROMWELL=/share/ClusterShare/software/contrib/micgea/cromwell/38/cromwell-38.jar  # requires java 8/1.8
CROMWELL=/share/ClusterShare/software/contrib/micgea/cromwell/68.1/cromwell-68.1.jar  # requires java 11
# Host, port, database and platform defaults
DBHOST="0.0.0.0"
DBPORT="40008"
DBNAME="cromwell"
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
                -M|--mysql)
                        MYSQL="$2"
                        shift
                        shift
                        ;;
                -R|--mysql_rundir)
                        MYSQL_RUNDIR="$2"
                        shift
                        shift
                        ;;
                -H|--dbhost)
                        DBHOST="$2"
                        shift
                        shift
                        ;;
                -P|--dbport)
                        DBPORT="$2"
                        shift
                        shift
                        ;;
                -n|--dbname)
                        DBNAME="$2"
                        shift
                        shift
                        ;;
                -c|--cromwell)
                        CROMWELL="$2"
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

echo "Database host       = ${DBHOST}"
echo "Database port       = ${DBPORT}"
echo "Database name       = ${DBNAME}"
echo "Cromwell port       = ${CROMPORT}"
echo "Platform            = ${PLATFORM}"
echo "Multi-sample mode   = ${MULTI}"
echo "MySQL directory     = ${MYSQL}"
echo "MySQL run directory = ${MYSQL_RUNDIR}"
echo "Cromwell location   = ${CROMWELL}"
if [[ -n $1 ]]; then
        echo "Remaining arguments: "
        echo "$@"
fi

if [ "${DRYRUN}" == "1" ]
then
        exit 0
fi

# Set cromwell basename variable
CROMWELL_BN="$(basename ${CROMWELL})"

# Configure mutect2.conf
# Set the MySQL hostname, port, and database name
sed -i -e "s/DBHOST_TO_SED/${DBHOST}/g" \
    -e "s/DBPORT_TO_SED/${DBPORT}/g" \
    -e "s/DBNAME_TO_SED/${DBNAME}/g" \
    -e "s/CROMWELL_PORT_TO_SED/${CROMPORT}/g" \
    ./workflow/mutect2.conf
# Set the GCP service account details
# Requires a service account email on the first line of .service_account.email.txt
# Requires a service account key supplied in .service_account.key.pem
if [ -f .service_account.email.txt ] && [ -f .service_account.key.pem ]
then
    EMAIL=$(awk 'NR==1' .service_account.email.txt | tr -d '\n')
    sed -i -e "s#ROOT_PATH_TO_SED#${PWD}#g" \
        -e "s/GCP_SA_EMAIL_TO_SED/${EMAIL}/g" \
        ./workflow/mutect2.conf
else
    echo "NO GCP SERVICE ACCOUNT DETAILS PROVIDED!"
    exit 0
fi

# Configure options files
# '#' delimiters are used to enable path substitution
sed -i -e "s#PWD_TO_SED#${PWD}#g" ./workflow/options.json
sed -i -e "s#PWD_TO_SED#${PWD}#g" -e "s/ID_TO_SED/${DBNAME}/g" ./workflow/options.google.json

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
sed -i -e "s/CROMWELL_PORT_TO_SED/${CROMPORT}/g" ./create_pon.sh
if [ "${PLATFORM}" == "GCP" ]
then
    sed -i -e "s/options\.json/options\.google\.json/g" ./run.sh
    sed -i -e "s/options\.json/options\.google\.json/g" ./run_chip_only.sh
    sed -i -e "s/options\.json/options\.google\.json/g" ./create_pon.sh
fi
sed -i -e "s/MULTI_TO_SED/${MULTI}/g" ./run.sh
sed -i -e "s/MULTI_TO_SED/${MULTI}/g" ./run_chip_only.sh
cd workflow
zip cromwell-kccg-mutect2.multi.dep.zip cromwell-kccg-mutect2.wdl cromwell-kccg-mutect2.chip.wdl
cd ..

# Configure start_cromwell.sh
ln -s ${CROMWELL}
sed -i -e "s/CROMWELL_JAR_TO_SED/${CROMWELL_BN}/g" ./start_cromwell.sh

# Create the database for Cromwell
echo 'CREATE DATABASE '"${DBNAME}"';
GRANT ALL PRIVILEGES ON '"${DBNAME}"'.* TO '\''cromwell_user'\''@'\''localhost'\'' WITH GRANT OPTION;
GRANT ALL PRIVILEGES ON '"${DBNAME}"'.* TO '\''cromwell_user'\''@'\''%'\'' WITH GRANT OPTION;' \
| ${MYSQL}/bin/mysql -u root --password=password --socket="${MYSQL_RUNDIR}"/socket
