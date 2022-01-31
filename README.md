# cromwell-kccg-mutect2

Run the KCCG GATK4 Mutect2 somatic variant calling pipeline using the Cromwell workflow engine.

## Usage:

This pipeline can be run in several different modes:
* Full pipeline: In this mode the entire pipeline will be run, starting with input BAM/CRAM files and performing somatic variant calling. Optionally, the VCF output can be further annotated with either Annovar or VEP. Additionally, CHIP variant detection can be run as a final stage of the pipeline, or can be run separately at a later time.
* Panel of Normals creation: This mode generates the Panel of Normals (PoN) from a set of input files. Input can either be BAM/CRAM alignment files, or pre-called VCF files.
* Annotation-only: The pipeline can be run in annotation-only mode to annotate pre-called VCF files. Both Annovar and VEP are available.
* CHIP-only: CHIP variant detection can be run stand-alone on pre-called VCF files.

The pipeline can also be run on several platforms:
* Local HPC cluster: Currently this supports the Sun Grid Engine (SGE) HPC cluster at the Garvan using the Cromwell workflow manager.
* Google Cloud Platform (GCP): This also uses the Cromwell workflow manager to manage running the pipeline remotely on GCP.
* Terra: Terra is a GCP-based front-end that uses Cromwell in the background to run workflows. This pipeline can be run natively on Terra.
* DNAnexus: DNAnexus is an Amazon Web Serivces (AWS)-based cloud platform. While it doesn't support WDL workflows natively, the dxCompiler software provided by DNAnexus can compile this script into a DNAnexus-ready workflow.

### Running on the Garvan HPC and Google Cloud Platform

For running the pipeline on the local cluster (SGE backend) or the Google Cloud Platform, a Cromwell server must be first set up on the local cluster. Cromwell handles orchestrating the pipeline and submitting each task to either the cluster or GCP.

Cromwell is not required for running the pipeline on DNAnexus or Terra. For instructions for those platforms, see further below.

#### Configuring a MySQL server
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

#### Configuring a Cromwell server

First, clone the repository to a suitable location.

```bash
# Clone the repository into your working directory
git clone https://git.gimr.garvan.org.au/micgea/cromwell-kccg-mutect2-chip.git
cd cromwell-kccg-mutect2
```

##### Using a Google Service Account

If running the pipeline on the Google Cloud Platform, before configuring the Cromwell server, Google Service Account credentials are required. Set up a Compute Engine Service Account and an authentication key (https://console.cloud.google.com/apis/credentials). Add the Service Account email to the first line in a new file within the top-level run directory with the name '.service_account.email.txt'. Similarly, add the key to a file with the name '.service_account.key.pem'. These strings can be found in the JSON file downloaded from GCP when setting up the authentication key.

```bash
echo "<PROJECT_ID>-compute@developer.gserviceaccount.com" > .service_account.email.txt
echo -e -n "-----BEGIN PRIVATE KEY-----\nABCDEFG<.....>TUVWXYZ\n-----END PRIVATE KEY-----\n" > .service_account.key.pem
chmod 600 .service_account.email.txt
chmod 600 .service_account.key.pem
```

Now configure the Cromwell server by running configure_cromwell.sh.

```bash
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
    -C  # Optional: enable call caching (off by default).
```

#### Running Cromwell server

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

#### Configuring workflow

Inside the input directory are several JSON files that describe the input files to the pipeline. One of these will be used for a given run, depending on which pipeline is being run:

| Filename | Description |
| -------- | ----------- |
| inputs.pon.json | Panel of normals creation |
| inputs.json | Full pipeline, single input mode |
| inputs.multi.json | Full pipeline, batch mode |
| inputs.chip.json | CHIP detection only, single input mode |
| inputs.chip.multi.json | CHIP detection only, batch mode |
| inputs.annovar.json | Annovar annotation only, single input mode |
| inputs.annovar.multi.json | Annovar annotation only, batch mode |
| inputs.vep.json | VEP annotation only, single input mode |
| inputs.vep.multi.json | VEP annotation only, batch mode |

Edit the appropriate JSON file in your favourite text editor. The input file is a series of key: value pairs. The templates provided have values of "REQUIRED_*" and "OPTIONAL_*" for required and optional fields, respectively, with with '*' indicating the type of input required (e.g. OPTIONAL_FILE or REQUIRED_STRING). If not using an optional parameter, simply delete the entire line. Some defaults are also provided and can be left as-is or changed if desired.

Now configure the run by running configure_run.sh.

```bash
# Run configuration script

# To see the available options and the defaults:
./configure_run.sh -h

# To perform a dry-run and only print the configuration settings without modifying any files:
./configure_run.sh -d [OPTIONS]

# When ready, run the configuration script with desired settings:
DB_NAME="example_run"  # Use only a-z, A-Z, and '_'. DO NOT use hyphen '-' character
CROMWELL_PORT=8007
PLATFORM=GCP

./configure_cromwell.sh \
    -n ${DB_NAME} \
    -p ${CROMWELL_PORT} \
    -f ${PLATFORM} \
    -m  # Optional: set the run to batch/multi-sample mode. This requries a TSV file describing all the input files.
```

##### A note about memory requirements

When running the workflow on the Garvan HPC, the requested memory will be multiplied by the requested number of cores. However, when running on GCP, Terra, or DNAnexus, the requested memory will be the total memory allocated to a job. To ensure the proper amount of memory is being requested in each case, set the "mem_per_core" parameter in the input JSON file, i.e. for the local HPC, set "mem_per_core" to true, and for GCP/Terra/DNAnexus, set it to false.

##### Prerequisite files

The following files are required to run parts of the pipeline.

| Parameter(s) | Description | Example |
| ------------ | ----------- | ------- |
| gnomad, gnomad_idx | Germline reference VCF and index containing common and rare variant population allele frequencies | https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz, https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi |
| variants_for_contamination, variants_for_contamination_idx | VCF and index containing common variants and allele frequencies for calculating contamination | https://storage.googleapis.com/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz, https://storage.googleapis.com/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz.tbi |
| vep_cache_archive | TAR.GZ archive file containing a VEP cache for annotating variants offline | http://ftp.ensembl.org/pub/release-103/variation/vep/homo_sapiens_vep_103_GRCh38.tar.gz |
| vep_loftee_ancestor_fa, vep_loftee_ancestor_fai, vep_loftee_ancestor_gzi | FASTA file and index for running VEP + LOFTEE | https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz, https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz.fai, https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz.gzi |
| vep_loftee_conservation_sql | PhyloCSF database for conservation filters in VEP + LOFTEE | https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/loftee.sql.gz |
| annovar_archive | TAR.GZ archive file containing the necessary files to run ANNOVAR | /share/ClusterShare/software/contrib/micgea/chip/data/annovar_files.tar.gz |
| whitelist_archive, whitelist_filter_archive | TAR.GZ archive file containing the necessary files to run the CHIP detection and variant whitelisting | /share/ClusterShare/software/contrib/micgea/chip/data/whitelist_filter_files.tar.gz |

##### Batch/multi-sample mode

If running the pipeline in batch/multi-sample mode, you will need to supply a file to the pipeline describing each sample. A template is provided in inputFile.tsv

For the full pipeline, the input file is a tab-separated values (TSV) file. Each line describes the files required to run one sample. The file must contain either TWO or FOUR tab-delimited columns. The columns describe the tumor BAM, tumor BAI, normal BAM (optional), and normal BAI (optional) files, respectively.

For the CHIP-only pipeline, the input file should contain two columns: the first containing the sample name (which must match the tumor sample name in the VCF header), and the file path to the filtered (and possibly annotated) VCF file to be used as input to the CHIP detection pipeline.

For annotating VCF files with VEP, a two-column input file is required: the file path to the VCF file, and the path to the VCF index.

For annotating with Annovar, only a single column is required: the file path to the VCF file.

#### Create a panel of normals (optional)

While it is possible to run the pipeline without a panel of normals (PoN), it is not recommended to do so. A generic PoN can be used; for example, the GATK best practices PoN generated from the 1000 genomes dataset: gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz. However, it is advisable to generate a PoN from a similar source to the samples being analysed, as this will help identify sequencing artefacts that may be called as false positive somatic variants. See https://gatk.broadinstitute.org/hc/en-us/articles/360035890631-Panel-of-Normals-PON- for further details on the best practices for generating a PoN. To generate a PoN, use the create_pon.sh script.

Similar to the batch/multi-sample mode, the panel of normals pipeline requires a two-column input TSV file describing the BAM and BAI file pairs for each sample. For reference, the panel of normals pipeline in effect runs the Mutect2 pipeline in tumor-only mode for each sample, hence the need for a two-column TSV file.

```bash
# Create an inputFile.tsv file. It should be in the following format:
cat inputFile.tsv

    /path/to/first/pon/sample.bam   /path/to/first/pon/sample.bai
    /path/to/second/pon/sample.bam   /path/to/secon/pon/sample.bai

# Ensure inputFile.tsv is specified in input/inputs.pon.json
cat input/inputs.pon.json

    {
        ...
        "Mutect2CHIP_Panel.bam_bai_list": "./inputFiles.tsv",
        ...
    }

# Create the PoN
./create_pon.sh
```

#### Running the workflow

Once everything is set up and configured, run the workflow as follows:

```bash
# Submit the workflow to Cromwell
# Run one of the following scripts:

# Option 1: run the full pipeline (BAM --> somatic variant calling --> CHIP detection)
./run.sh

# Option 2: only run the CHIP detection stage
./run_chip_only.sh

# Option 3: only run VEP or Annovar annotation
./run_vep_only.sh
./run_annovar_only.sh
```

You can check on the state of the run by re-attaching the Cromwell screen session:

```bash
screen -r "CromwellPort${CROMWELL_PORT}"
```

#### Aborting the workflow

To abort the workflow, simply run the following:

```bash
python3 abort.py
```

Alternatively, you can directly run the abort script by passing it the run ID. This is the ID returned in JSON format when submitting the workflow. It can be found in run_id.txt.

```bash
cat run_id.txt

    {"id":"3deb2054-445f-400a-afcd-465a50fc44ba","status":"Submitted"}

./abort.sh "3deb2054-445f-400a-afcd-465a50fc44ba"
```

#### Resetting the database

If you need to start from scratch, you can choose to delete and re-create the Cromwell database by running reset_database.sh. This may be useful if the workflow fails to abort properly.

```bash
./reset_database.sh
```

### Running on Terra

#### TODO: finish this section

Running the pipeline on Terra is similar to running on GCP, except that a Cromwell server is not required. Instead, all that is required is to:

1. Upload the relevant workflow up to the Terra platform.
2. Generate the appropriate input JSON file as described above.
    1. Ensure that all files are uploaded to a GCP bucket and that their gs:// URLs are entered into the input JSON file.
3. Upload the input JSON file to Terra.
5. Run the workflow on Terra.

### Running on DNAnexus

Running the pipeline on DNAnexus (and DNAnexus-based platforms such as the UK Biobank Reasearch Analysis Platform (RAP)) consists of four steps: compiling the pipeline up to DNAnexus, uploading data to a project, configuring the pipeline, and running the pipeline.

Ensure that you have installed [dx-toolkit](https://documentation.dnanexus.com/downloads) and that the ```dx``` command is available in your path.

#### Compiling the pipeline

Compilation of the workflow requires the dxCompiler tool provided by DNAnexus: [dxCompiler](https://github.com/dnanexus/dxCompiler). The latest version ([v2.8.3](https://github.com/dnanexus/dxCompiler/releases/tag/2.8.3)) requires Java v11 to be installed.

Prior to compilation, ensure that you are logged in to your DNAnexus account, have created a project for the analysis, and have selected that project with ```dx```.

```bash
# If you already have a project ready to go, dx login will let you interactively select it from a list
dx login

# If you first need to create a project, run these commands instead
dx login --noprojects
dx new project --name <PROJECT NAME> --region <REGION> --bill-to <BILLING ACCOUNT>
dx select --name <PROJECT NAME>
```

Once ready, change into the ```workflow/``` directory and compile the workflow as follows:
```bash
cd workflow/
java -jar /PATH/TO/dxCompiler-2.8.3.jar compile cromwell-kccg-mutect2.wdl -destination <PROJECT NAME>:/PATH/TO/PLACE/WORKFLOW/ -f -reorg
# The -f option will force dxCompiler to overwrite the workflow and/or applets if they are already present in the workflow directory on DNAnexus.
# The -reorg option will separate the intermediate files away from the main workflow outputs
```

If this completes successfully, the workflow will be ready to run, and located in the selected project at the location defined with the ```-destination``` parameter.

#### Uploading data to DNAnexus

Upload all necessary data to the DNAnexus project using dx-toolkit:

```bash
dx upload --path <PROJECT NAME>:/PATH/TO/UPLOAD/DATA [-r] filename [filename ...]

# -r means upload a directory recursively, similar to cp -r
```

For example, to upload a single file named ```test.txt``` to ```test_project``` at the location ```/data/```:

```bash
dx upload --path test_project:/data/ test.txt
```

Or to upload an entire directory called ```test_dir/```:

```bash
dx upload --path -r test_project:/data/ test_dir
```

#### Configuring the pipeline for DNAnexus

Use the input templates in the ```input/``` directory to configure the pipeline, similar to what is described above for Cromwell workflows. The only difference is to provide the DNAnexus paths for files. For example:

```json
{
    "Mutect2CHIP.intervals": "project-G7Q2GLb5FXyKGpYj4SpgKfk3:/data/exome_evaluation_regions.v1.interval_list"
    ...
}
```

#### Running the pipeline on DNAnexus

Finally, to run the pipeline, first ensure you have Python 3 installed and available on your PATH. Then, run the following:

```bash
./run_dx.sh -i input/<INPUT JSON FILE> -o <PROJECT NAME>:/PATH/TO/ANALYSIS/OUTPUT -w <PROJECT NAME>:/PATH/TO/WORKFLOW
```

For example, to run the main, single-sample pipeline loated in the DNAnexus project ```project-G7Q2GLb5FXyKGpYj4SpgKfk3``` at ```/workflows/single_sample/Mutect2CHIP``` using the configuration in ```input/inputs.json```, and to output the data to ```/output/run_1```:

```bash
./run_dx.sh -i input/inputs.json -o project-G7Q2GLb5FXyKGpYj4SpgKfk3:/output/run_1/ -w project-G7Q2GLb5FXyKGpYj4SpgKfk3:/workflows/single_sample/Mutect2CHIP
```
