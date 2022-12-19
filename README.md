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

#### Setting up a Cromwell Server

We currently have a production Cromwell server in development, the repository for which can be found [here](https://git.gimr.garvan.org.au/CCG/kccg-cromwell/-/tree/dev). Follow the steps in the README for that repository to set up a Cromwell server for running this pipeline.

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
RUN_NAME="example_run"  # Use only a-z, A-Z, and '_'. DO NOT use hyphen '-' character
CROMWELL_PORT=8007
PLATFORM=GCP

./configure_run.sh \
    -n ${RUN_NAME} \
    -p ${CROMWELL_PORT} \
    -f ${PLATFORM} \
    -m  # Optional: set the run to batch/multi-sample mode. This requries a TSV file describing all the input files.
```

##### Prerequisite files

The following files are required to run parts of the pipeline.

| Parameter(s) | Description | Example |
| ------------ | ----------- | ------- |
| gnomad, gnomad_idx | Germline reference VCF and index containing common and rare variant population allele frequencies | https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz, https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi |
| variants_for_contamination, variants_for_contamination_idx | VCF and index containing common variants and allele frequencies for calculating contamination | https://storage.googleapis.com/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz, https://storage.googleapis.com/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz.tbi |
| vep_cache_archive | TAR.GZ archive file containing a VEP cache for annotating variants offline | http://ftp.ensembl.org/pub/release-103/variation/vep/homo_sapiens_vep_103_GRCh38.tar.gz |
| vep_loftee_ancestor_fa, vep_loftee_ancestor_fai, vep_loftee_ancestor_gzi | FASTA file and index for running VEP + LOFTEE | https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz, https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz.fai, https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz.gzi |
| vep_loftee_conservation_sql | PhyloCSF database for conservation filters in VEP + LOFTEE | https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/loftee.sql.gz |
| annovar_archive | TAR.GZ archive file containing the necessary files to run ANNOVAR | https://storage.cloud.google.com/kccg-somvar-data/annovar_files.tar.gz |
| whitelist_archive, whitelist_filter_archive | TAR.GZ archive file containing the necessary files to run the CHIP detection and variant whitelisting | https://storage.cloud.google.com/kccg-somvar-data/whitelist_filter_files.tar.gz |

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
./run.sh -m pon
```

#### Running the workflow

Once everything is set up and configured, run the workflow as follows:

```bash
# Submit the workflow to Cromwell
# Run the pipeline in one of the following modes:

# Option 1: run the full pipeline (BAM --> somatic variant calling --> CHIP detection)
RUNMODE=full

# Option 2: only run the CHIP detection stage
RUNMODE=chip

# Option 3: only run VEP or Annovar annotation
RUNMODE=annovar
# RUNMODE=vep

# Then execute the run script in the appropriate mode
./run.sh -m ${RUNMODE}
```

#### Workflow status and termination

To check on the workflow's status, run the script `status.sh`. Similarly, to terminate the workflow, run `abort.sh`. These scripts both pull the run ID out of run_id.txt, which is generated upon job submission with run.sh. Alternatively, these scripts can be provided the run ID directly as the first argument.

```bash
# Check the run ID
cat run_id.txt

    {"id":"3deb2054-445f-400a-afcd-465a50fc44ba","status":"Submitted"}

# Get workflow status
./status.sh  # Alternatively, ./status.sh 3deb2054-445f-400a-afcd-465a50fc44ba

    {"status":"Running","id":"3deb2054-445f-400a-afcd-465a50fc44ba"}

# Terminate the workflow
./abort.sh  # Alternatively, ./abort.sh 3deb2054-445f-400a-afcd-465a50fc44ba

    {"id":"3deb2054-445f-400a-afcd-465a50fc44ba","status":"Aborting"}
```

### Running on Terra

Running the pipeline on Terra is similar to running on GCP, except that a Cromwell server is not required. All that is required is to:

1. Create a workspace on Terra
2. Add the workflow(s) to the Terra workspace.
3. Upload all requisite files to a GCP bucket.
    1. If running in batch mode, ensure the inputFiles.tsv files is uploaded to GCP as well.
    2. Ensure all input files are recorded in the input JSON file with their gs:// URLs.
4. Upload the input JSON file to the Terra workspace.
5. Run the workflow on Terra.

#### Terra workspace

This workflow is already hosted on Terra in a [private workspace](https://app.terra.bio/#workspaces/terra-kccg-production/terra-kccg-somvar-pipeline). If you have access to this workspace, you can clone it and start running the pipeline (see below).

For setting up your own workspace, see [this Google Docs document](https://docs.google.com/document/d/1pZVTxjRJfAyWYiFWmmF_o31lORQ8a-xmqdK1zE3RDyk) for detailed instructions on setting up a Terra workspace. In short, a billing account on GCP is required to be linked with a billing account on Terra, to which the new workspace can be added. A Terra group must also be created, containing any Terra users that wish to run the pipeline, and the group's firecloud email added as a "Storage Object Creator" and "Storage Object Viewer" on any GCP buckets used to host data as part of the workflow.

Once set up, the workflows can be added to Terra via Dockstore. Dockstore monitors a public git repository containing the workflows and supplies these workflows for Terra to use. The public repository of this pipeline is available [here](https://github.com/mgeaghan-garvan/cromwell-kccg-mutect2-chip), and is already hosted on Dockstore for use with Terra. Note that this repository will typically only be updated with the main, production-ready branch, and may not contain updates in the development branch. To add the workflow to Terra:

1. Go to the Terra workspace.
2. Go to the "Workflows" tab.
3. Click on "Find a Workflow".
4. Click on "Dockstore" under "Find Additional Workflows".
5. Search for the workflows with "github.com/mgeaghan-garvan/cromwell-kccg-mutect2-chip". Click on the link to the desired workflow.
6. Click on "Terra" under the "Launch with" section.
7. Set the workflow name and pick a destination workspace, then click "Import".

#### Setting up and running a workflow

As mentioned above, first ensure that all requisite files, including the inputFiles.tsv containing paths to all BAMs/CRAMs/VCFs for batch runs are present in a GCP bucket to which the Terra workspace has access. Once data has been uploaded to GCP:

1. Go to the Terra workspace.
2. Go to the "Workflows" tab.
3. Click on the relevant workflow.
4. Under the "Inputs" tab, drag and drop the input JSON file, or search for it by clicking on "upload json". Alternatively, the input fields can be filled out manually.
5. Click on "Run Analysis".

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

Ensure you have Python 3 installed and available on your PATH. Then, run the configuration script to create the DNAnexus submission script(s).

```bash
./configure_run.sh \
    -f DX \
    -w <DNAnexus workflow path> \
    -i <Input JSON file> \
    -o <DNAnexus output path> \
    -b <Local input TSV file>  # This is optional (see below)
```

##### Submitting multiple samples to DNAnexus
You can use the ```-b``` flag for configure_run.sh to configure the run to submit each sample in a batch as a separate job to DNAnexus. In this way, if a sample fails, the other samples won't be affected. The argument to the -b flag is the same inputFiles.tsv file you would supply for a batch run on Cromwell. NOTE: when running the pipeline in this way, you need to supply the corresponding single sample JSON file and DNAnexus workflow to the ```-i``` flag, since the pipeline is actually running multiple single-sample mode runs in parallel. For example:

```bash
./configure_run.sh \
    -f DX \
    -w project-X7Gjf3BJbVIJABcd1234KI95:/workflow/Mutect2CHIP \
    -i ./input/inputs.json \
    -o project-X7Gjf3BJbVIJABcd1234KI95:/output \
    -b ./input/inputFiles.tsv
```

#### Running the pipeline on DNAnexus

Finally, to run the pipeline, run the following:

```bash
./run.sh -m dx
```
