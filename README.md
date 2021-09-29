# cromwell-kccg-mutect2

Run the KCCG GATK4 Mutect2 somatic variant calling pipeline using the Cromwell workflow engine.

## Usage:

### Configuring MySQL server
If there is no MySQL server currently running on the cluster, run the following commands:
```bash
# Choose an empty port number.
MYSQL=/home/glsai/mysql/mysql-5.7.27-linux-glibc2.12-x86_64/
MYSQL_RUNDIR=/home/glsai/mysql/mysql-5.7.27-linux-glibc2.12-x86_64/
sqlport=40008
 
screen -S MySQLPort${sqlport}
 
# Unfortunately have to run these lines again inside screen session
MYSQL=/home/glsai/mysql/mysql-5.7.27-linux-glibc2.12-x86_64/
MYSQL_RUNDIR=/home/glsai/mysql/mysql-5.7.27-linux-glibc2.12-x86_64/
sqlport=40008
 
${MYSQL}/bin/mysqld --no-defaults --user=glsai --datadir=$MYSQL_RUNDIR/data/ --basedir=$MYSQL_RUNDIR/ --log-error=$MYSQL_RUNDIR/log/mysql.err --pid-file=$MYSQL_RUNDIR/mysql.pid --socket=$MYSQL_RUNDIR/socket --port=${sqlport}
 
# Then exit screen session once our SQL is running. When its running it says or does nothing.
# Ctrl + A and D
 
${MYSQL}/bin/mysql -u root -p --socket=$MYSQL_RUNDIR/socket
# password: password
 
# if cromwell_user not already created else skip these two lines
# you can check if it already exists with the following query: "SELECT user FROM mysql.user;"
CREATE USER 'cromwell_user'@'localhost' IDENTIFIED BY '123456781!aA';
CREATE USER 'cromwell_user'@'%' IDENTIFIED BY '12345678!1aA';
```

### Configuring Cromwell server

```bash
# Clone the repository into your working directory
git pull https://git.gimr.garvan.org.au/micgea/cromwell-kccg-mutect2.git
cd cromwell-kccg-mutect2

# Run configuration script

# To see the available options and the defaults:
./configure_cromwell.sh -h

# To perform a dry-run and only print the configuration settings without modifying any files:
./configure_cromwell.sh -d [OPTIONS]

# When ready, run the configuration script with desired settings:
DB_HOST="0.0.0.0"  # Default. Needs to be '0.0.0.0' rather than 'localhost'
DB_PORT=40008
DB_NAME="example_run"  # Use only a-z, A-Z, and '_'. DO NOT use hyphen '-' character
CROMWELL_PORT=8007
PLATFORM=GCP
MYSQL_DIR=/home/glsai/mysql/mysql-5.7.27-linux-glibc2.12-x86_64/
MYSQL_RUN_DIR=/home/glsai/mysql/mysql-5.7.27-linux-glibc2.12-x86_64/
CROMWELL_JAR=/share/ClusterShare/software/contrib/micgea/cromwell/68.1/cromwell-68.1.jar  # Cromwell JAR file location

./configure_cromwell.sh \
    -H ${MYSQL_HOST} \
    -P ${DB_PORT} \
    -n ${DB_NAME} \
    -p ${CROMWELL_PORT} \
    -f ${PLATFORM} \
    -M ${MYSQL_DIR} \
    -R ${MYSQL_RUN_DIR} \
    -c ${CROMWELL_JAR} \
    -m  # Optional: set the run to batch/multi-sample mode. This requries a TSV file describing all the input files.
```

### Running Cromwell server

```bash
# Create a screen session
screen -S "CromwellPort${CROMWELL_PORT}"

# Ensure that Java 11 is the active Java version
java -version
# If it isn't check the currently active modules
module list
# If an older version of Java is loaded, unload it. This should return the default Java version to 8/v1.8
module unload centos6.10/ccg/java/1.7.0_25
# Now load Java 11
module load centos6.10/shacar/java/jdk-11.0.2

# Start the Cromwell server
./start_cromwell.sh
```

Once running, you can detatch the session with Ctrl + A, then D.

### Configuring workflow

Inside the workflow directory are several JSON files that describe the input files to the pipeline. One of these will be used for a given run, depending on which pipeline is being run:

    Panel of normals creation:                     inputs.pon.json
    Full pipeline, single input mode:              inputs.json
    Full pipeline, batch/mulit-sample mode:        inputs.multi.json
    CHIP detection-only, single input mode:        inputs.chip.json
    CHIP detection-only, batch/multi-sample mode:  inputs.chip.multi.json

Edit the appropriate JSON file in your favourite text editor. The input file is a series of key: value pairs. The templates provided have values of "REQUIRED" and "OPTIONAL" for required and optional fields, respectively. If not using an optional parameter, simply delete the entire line. Some defaults are also provided and can be left as-is or changed if desired.
 
#### Batch/multi-sample mode

If running the pipeline in batch/multi-sample mode, you will need to supply a file to the pipeline describing each sample. A template is provided in inputFile.tsv

For the full pipeline, the input file is a tab-separated values (TSV) file. Each line describes the files required to run one sample. The file must contain either TWO or FOUR tab-delimited columns. The columns describe the tumor BAM, tumor BAI, normal BAM (optional), and normal BAI (optional) files, respectively.

For the CHIP-only pipeline, the input file should contain only a single column, with each line describing the filtered (and possibly annotated) VCF file for a single sample to be used as input to the CHIP detection pipeline.

### Create a panel of normals (optional)

If a panel of normals is required, run the create_pon.sh script.

Similar to the batch/multi-sample mode, the panel of normals pipeline requires a two-column input TSV file describing the BAM and BAI file pairs for each sample. For reference, the panel of normals pipeline in effect runs the Mutect2 pipeline in tumor-only mode for each sample, hence the need for a two-column TSV file.

```bash
# Create an inputFile.tsv file. It should be in the following format:
cat inputFile.tsv

    /path/to/first/pon/sample.bam   /path/to/first/pon/sample.bai
    /path/to/second/pon/sample.bam   /path/to/secon/pon/sample.bai

# Ensure inputFile.tsv is specified in workflow/inputs.pon.json
cat workflow/inputs.pon.json

    {
        ...
        "Mutect2CHIP_Panel.bam_bai_list": "./inputFiles.tsv",
        ...
    }

# Create the PoN
./create_pon.sh
```

### Running the workflow

Once everything is set up and configured, run the workflow as follows:

```bash
# Submit the workflow to Cromwell
# Run one of the following scripts:

# Option 1: run the full pipeline (BAM --> somatic variant calling --> CHIP detection)
./run.sh

# Option2: only run the CHIP detection stage
./run_chip_only.sh
```

You can check on the state of the run by re-attaching the Cromwell screen session:

```bash
screen -r "CromwellPort${CROMWELL_PORT}"
```
