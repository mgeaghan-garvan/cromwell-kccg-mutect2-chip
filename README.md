# cromwell-kccg-mutect2

Run the KCCG GATK4 Mutect2 somatic variant calling pipeline using the Cromwell workflow engine.

## Usage:

This pipeline can be run in several different modes:
* Full pipeline: In this mode the entire pipeline will be run, starting with input BAM/CRAM files and performing somatic variant calling with GATK Mutect2. Optionally, the VCF output can be further annotated with Annovar, VEP, and/or SpliceAI. Finally, CHIP variant annotation can be run as a final stage of the pipeline.
* Panel of Normals creation: This mode generates the Panel of Normals (PoN) from a set of input files. Input can either be BAM/CRAM alignment files, or pre-called VCF files.
* Annotation-only: The pipeline can be run in annotation-only mode to annotate pre-called VCF files. Annovar, VEP, and SpliceAI are available.
* CHIP-only: CHIP variant detection can be run stand-alone on pre-called VCF files.

A cohort-wide CHIP filter workflow is also available, which is recommended to be run following CHIP annotation. This cohort-wide stage detects and filters out mutations that appear above a specified threshold in the cohort, as well as mutations that are multiallelic across the cohort.

This repository also includes a cohort annotation workflow, which can significantly reduce the costs of annotating large cohorts by assembling a single sites-only VCF that can be annotated with Annovar, VEP, or SpliceAI; these annotations are then added to the original VCFs using bcftools, thus reducing the potential overhead of annotating the same variant multiple times in different samples.

The pipeline can also be run on several platforms:
* Local HPC cluster: Currently this supports the Sun Grid Engine (SGE) HPC cluster at the Garvan using the Cromwell workflow manager.
* Google Cloud Platform (GCP): This also uses the Cromwell workflow manager to manage running the pipeline remotely on GCP.
  * The pipeline can also be run on on the Cromwell server hosted by the Centre for Population Genomics (CPG). This involves submitting the workflow to the CPG's hail batch server via their [analysis-runner](https://github.com/populationgenomics/analysis-runner) tool.
* Terra: Terra is a GCP-based front-end that uses Cromwell in the background to run workflows. This pipeline can be run natively on Terra.
* DNAnexus: DNAnexus is an Amazon Web Serivces (AWS)-based cloud platform. While it doesn't support WDL workflows natively, the dxCompiler software provided by DNAnexus can compile this script into a DNAnexus-ready workflow.
  * A single-job version of the Mutect2 workflow (no annotation or CHIP detection) has been created specifically for running on DNAnexus, as it runs more efficiently on that platform compared to multiple smaller jobs.

### Running on the Garvan HPC and Google Cloud Platform

For running the pipeline on the local cluster (SGE backend) or the Google Cloud Platform, a Cromwell server must be first set up on the local cluster. Cromwell handles orchestrating the pipeline and submitting each task to either the cluster or GCP.

Cromwell is not required for running the pipeline on DNAnexus or Terra. For instructions for those platforms, see further below.

#### Setting up a Cromwell Server

We currently have a production Cromwell server in development, the repository for which can be found [here](https://git.gimr.garvan.org.au/CCG/kccg-cromwell/-/tree/dev). Follow the steps in the README for that repository to set up a Cromwell server for running this pipeline.

#### Configuring a workflow

Inside the input directory are several JSON files that describe the input files to the pipeline. One of these will be used for a given run, depending on which pipeline is being run:

| Filename | Description |
| -------- | ----------- |
| inputs.pon.json | Panel of normals creation |
| inputs.m2.json | Mutect2-only pipeline |
| inputs.m2.sj.json | Mutect2-only pipeline, single-job version (for running on DNAnexus) |
| inputs.chip.json | CHIP detection only |
| inputs.chip_cohort.json | CHIP cohort-wide filter |
| inputs.annovar.json | Annovar annotation only |
| inputs.vep.json | VEP annotation only |
| inputs.spliceai.json | SpliceAI annotation only |
| inputs.full.json | Full pipeline, including Mutect2 and optionally one or more of Annovar, VEP, SpliceAI, and CHIP annotation |
| inputs.annotate_cohort.json | Annovar/VEP/SpliceAI annotation of an entire cohort |

Edit the appropriate JSON file in your favourite text editor. The input file is a series of key: value pairs. The templates provided have values of "REQUIRED_*" and "OPTIONAL_*" for required and optional fields, respectively, with with '*' indicating the type of input required (e.g. OPTIONAL_FILE or REQUIRED_STRING). If not using an optional parameter, simply delete the entire line. Some defaults are also provided and can be left as-is or changed if desired.

See the tables in the [Input parameters](#input-parameters) section below for details on the parameters for each pipeline mode.

##### Workflow options

A second JSON file is required to instruct Cromwell how to run the pipeline, e.g. where to store the outputs, and whether to keep intermediate files after the run has finished. Two example JSON files have been provided which can be used with [the local SGE cluster](workflow_options/options.sge.json) and [Google Cloud](workflow_options/options.gcp.json). Placeholder values have been provided in these files in uppercase.

For further details on the parameters to this file, see [the Cromwell documentation](https://cromwell.readthedocs.io).

#### Running a workflow

A handy run script called [run.sh](run.sh) has been written to simplify submitting and monitoring the workflow when submitted directly to a Cromwell server. It can also be used to submit jobs to DNAnexus, which is described further down. The list of parameters to this script can be seen by running `./run.sh -h`:

```bash
./run.sh -h

Usage: ./manage.sh [MODE] [--id WORKFLOW_ID] [--host CROMWELL_HOST] [--port CROMWELL_PORT] [--wdl WDL_PATH] [--input INPUT_JSON] [--batch BATCH_INPUT_TSV] [--options OPTIONS_JSON] [--dxworkflow DX_WORKFLOW] [--dxoutput DX_OUTPUT] [--dxpriority DX_PRIORITY]
Display this help message: ./run.sh -h

        MODE:                       Either 'submit', 'dxsubmit', 'status', or 'abort'.           REQUIRED
        --id WORKFLOW_ID:           Workflow ID for which to get status or abort.                REQUIRED for 'status' and 'abort' modes.
        --host CROMWELL_HOST:       Host on which Cromwell web service is running.               Default: 'localhost'
        --port CROMWELL_PORT:       Port on which Cromwell web service is running.               Default: '8007'
        --wdl WDL_PATH:             Path to WDL file describing workflow to run.                 REQUIRED for 'submit' mode
        --input INPUT_JSON:         Path to JSON file describing workflow inputs.                REQUIRED for 'submit' and 'dxsubmit' mode
        --batch BATCH_INPUT_TSV:    Path to TSV file describing inputs for batch runs.           REQUIRED for 'dxsubmit' mode; OPTIONAL for 'submit' mode; if supplied, a batch job will be submitted
        --options OPTIONS_JSON:     Path to JSON file describing workflow options.               REQUIRED for 'submit' mode
        --dxworkflow DX_WORKFLOW:   DNAnexus path to workflow.                                   REQUIRED for 'dxsubmit' mode
        --dxoutput DX_OUTPUT:       DNAnexus path to place output.                               REQUIRED for 'dxsubmit' mode
        --dxpriority DX_PRIORITY:   DNAnexus run priority; either 'low', 'normal', or 'high'.    Default: 'low'
```

To submit a job to a Cromwell server, you need the host address and the port on which the server is running, and the paths to the WDL workflow, input JSON file, and workflow options JSON file. Submission can be performed with:

```bash
# For example, to run the full pipeline
./run.sh submit --host <CROMWELL_HOST> --port <CROMWELL_PORT> --wdl workflow/full.wdl --input input/config/inputs.full.json --options workflow_options/options.sge.json
```

A successful workflow submission will return a workflow ID which can be used to monitor and terminate a running workflow:

```bash
# Check up on a run
./run.sh status --host <CROMWELL_HOST> --port <CROMWELL_PORT> --id <WORKFLOW_ID>

# Terminate a run
./run.sh abort --host <CROMWELL_HOST> --port <CROMWELL_PORT> --id <WORKFLOW_ID>
```

#### Running a batch of samples

A batch of jobs can be submitted all at once to a Cromwell server by supplying a TSV file to the `--batch` parameter of `run.sh`. The structure of the TSV file differs slightly between the different workflows, as different workflows have different input file requirements. The workflows that support a batch mode are:

- The full pipeline: [full.wdl](workflow/full.wdl)
- Mutect2: [m2.wdl](workflow/m2.wdl)
- All annotation-only pipelines: [annovar.wdl](workflow/annovar.wdl), [vep.wdl](workflow/vep.wdl) and [spliceai.wdl](workflow/spliceai.wdl)
- The CHIP-only pipeline: [chip.wdl](workflow/chip.wdl)

Template TSV files are supplied in `input/inputFiles/batch`:

- Full pipeline and Mutect2 template: [inputFiles.tsv](input/inputFiles/batch/inputFiles.tsv)
- Annovar and SpliceAI workflow template: [inputFiles.annotate.tsv](input/inputFiles/batch/inputFiles.annotate.tsv)
- VEP workflow template: [inputFiles.vep.tsv](input/inputFiles/batch/inputFiles.vep.tsv)
- CHIP-only workflow template: [inputFiles.chip.tsv](input/inputFiles/batch/inputFiles.chip.tsv)

For example, to run a batch of samples through the full pipeline, you would first create an input TSV file like so:

```bash
cat input/inputFiles/inputFiles.tsv

    sample_id   tumor_reads   tumor_reads_index   normal_reads  normal_reads_index
    SAMPLE_1    TUMOR_SAMPLE_1.bam    TUMOR_SAMPLE_1.bam.bai    NORMAL_SAMPLE_1.bam    NORMAL_SAMPLE_1.bam.bai
    SAMPLE_2    TUMOR_SAMPLE_2.bam    TUMOR_SAMPLE_2.bam.bai    NORMAL_SAMPLE_2.bam    NORMAL_SAMPLE_2.bam.bai
    TUMOR_ONLY_SAMPLE    TUMOR_ONLY_SAMPLE.bam TUMOR_ONLY_SAMPLE.bam.bai
```

Note that the header line contains the exact names of the workflow parameters. It also includes an additional first column for the sample ID. Also, for the full pipeline and Mutect2-only pipelines, the final two columns for normal reads are optional, as these parameters are optional in the pipeline.

Next, you would then submit the batch job like so:

```bash
./run.sh submit --host <CROMWELL_HOST> --port <CROMWELL_PORT> --wdl workflow/full.wdl --input input/config/inputs.full.json --batch input/inputFiles/inputFiles.tsv --options workflow_options/options.sge.json
```

The run script will create a copy of the input JSON file for each sample and automatically replace each sample-specific parameter with the corresponding value from the matching column.

The PoN generation workflow, cohort annotation workflow, and cohort-wide CHIP filter workflow are already native batch workflows that take a TSV file as an argument to one of their input parameters. As such, these do not work with the `--batch` mode of `run.sh`; instead, these TSV files should be passed to the corresponding parameter in their respective JSON configuration files. Template input TSV files for these workflows are also provided:

- Cohort annotation workflow: [inputFiles.annotate_cohort.tsv](input/inputFiles/batch/inputFiles.annotate_cohort.tsv)
- Cohort-wide CHIP filter workflow: [inputFiles.chip_cohort.tsv](input/inputFiles/batch/inputFiles.chip_cohort.tsv)
- PoN generation workflow - see [Creating a panel of normals (optional)](#creating-a-panel-of-normals-optional) below

#### Creating a panel of normals (optional)

While it is possible to run the pipeline without a panel of normals (PoN), it is not recommended to do so. A generic PoN can be used; for example, the GATK best practices PoN generated from the 1000 genomes dataset: gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz. However, it is advisable to generate a PoN from a similar source to the samples being analysed, as this will help identify sequencing artefacts that may be called as false positive somatic variants. See https://gatk.broadinstitute.org/hc/en-us/articles/360035890631-Panel-of-Normals-PON- for further details on the best practices for generating a PoN. To generate a PoN, use the create_pon.sh script.

The panel of normals pipeline requires a two-column input TSV file describing the BAM and BAI file pairs for each sample. It is also possible to supply a list of gzipped VCFs and their TBI index files, to run the PoN creation directly from VCFs. At least one of these files is necessary, and both can be provided in the case where some PoN-ready VCFs have already been generated.

Example TSV files are provided at [input/inputFiles/pon/inputFiles.pon.bams.tsv](input/inputFiles/pon/inputFiles.pon.bams.tsv) and [input/inputFiles/pon/inputFiles.pon.vcfs.tsv](input/inputFiles/pon/inputFiles.pon.vcfs.tsv). If a BAM TSV file is being provided, it needs to be passed to the `Mutect2_Panel.bam_bai_list` parameter in the input JSON file; if a VCF TSV file is being provided, it needs to be passed to `Mutect2_Panel.vcf_idx_list`.

### Running on the Centre for Population Genomics (CPG) Cromwell server

These workflows can also be run through the Cromwell server hosted by the CPG. This is required for processing datasets under their control. This will require the following:

- Installing the CPG's analysis-runner tool: [https://github.com/populationgenomics/analysis-runner](https://github.com/populationgenomics/analysis-runner)
- Cloning the CPG's fork of this pipeline: [https://github.com/populationgenomics/cromwell-kccg-mutect2-chip](https://github.com/populationgenomics/cromwell-kccg-mutect2-chip)

It is important to use `main` branch of the CPG's fork of this repository, as their Google Cloud infrastructure is set up to only run code that has been reviewed, approved, and merged into a main branch of a repository under their control.

The `analysis-runner` tool submits jobs to the CPG Hail Batch server. A Hail Batch job can be created to in turn submit a WDL pipeline to their Cromwell server. Currently, we have Hail Batch scripts set up under `scripts/cpg_cromwell` to submit the following workflows:

- Full analysis pipeline: [full.submit.py](scripts/cpg_cromwell/full.submit.py)
- Mutect2-only pipeline: [m2.submit.py](scripts/cpg_cromwell/m2.submit.py)
- PoN generation: [pon.submit.py](scripts/cpg_cromwell/pon.submit.py)

Since these workflows need to be submitted through Hail Batch, the configuration is slightly different, requiring a TOML file, rather than a JSON file. Additionally, batch submissions are not yet supported by these Hail Batch scripts. Template TOML files for the supported workflows are provided in the `input/config` directory. Each TOML file must start with a `workflow` section, describing the dataset, access level, and name of the workflow. More information about the workflow access levels can be found in the [CPG's documentation](https://github.com/populationgenomics/team-docs/tree/main/storage_policies#analysis-runner). The parameters to be passed to the WDL pipeline need to be defined in a separate `mutect2_chip` section. For example, the TOML file for a full workflow run will look something like the following:

```toml
[workflow]
dataset = 'kccg-genomics-med'
access_level = 'full'
name = 'workflow_name'

[mutect2_chip]
intervals = 'gs://path/to/intervals.interval_list'
ref_fasta = 'gs://path/to/hg38.fasta'
ref_fai = 'gs://path/to/hg38.fasta.fai'
ref_dict = 'gs://path/to/hg38.dict'
tumor_reads = 'gs://path/to/sample.cram'
tumor_reads_index = 'gs://path/to/sample.cram.crai'
normal_reads = 'gs://path/to/normal.cram'
normal_reads_index = 'gs://path/to/normal.cram.crai'
pon = 'gs://path/to/pon.vcf.gz'
pon_idx = 'gs://path/to/pon.vcf.gz.tbi'
...
```

To submit a workflow, run `analysis-runner` like so:

```bash
analysis-runner \
  --dataset kccg-genomics-med \
  --description "Example workflow" \
  --output-dir "relative/path/for/outputs" \
  --access-level full \
  --config input/config/inputs.full.toml \
  --image australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_workflows:latest \
  python3 scripts/cpg_cromwell/full.submit.py
```

The Cromwell server will write all final outputs to `gs://cpg-<DATASET>-main/mutect2-chip/<WORKFLOW_NAME>/`, where `<DATASET>` and `<WORKFLOW_NAME>` correspond to the strings given to the `dataset` and `name` parameters in the `workflow` section of the TOML file.

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

As mentioned above, first ensure that all requisite files, including the inputFiles.tsv containing paths to all BAMs/CRAMs/VCFs for batch runs (e.g. PoN generation, cohort-wide CHIP filter, cohort annotation) are present in a GCP bucket to which the Terra workspace has access. Once data has been uploaded to GCP:

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
java -jar /PATH/TO/dxCompiler-2.11.0.jar compile cromwell-kccg-mutect2.wdl -destination <PROJECT NAME>:/PATH/TO/PLACE/WORKFLOW/ -f -reorg
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

#### Submitting the pipeline to DNAnexus

Use the input JSON templates in the `input/config` directory to configure the pipeline, similar to what is described above for Cromwell workflows. The only difference is to provide the DNAnexus paths for files. For example:

```json
{
    "Mutect2CHIP.intervals": "project-G7Q2GLb5FXyKGpYj4SpgKfk3:/data/exome_evaluation_regions.v1.interval_list"
    ...
}
```

Similar to Cromwell runs, DNAnexus runs can be run either as a single job or in batch mode with the `--batch` parameter and an input TSV file. To run a workflow on DNAnexus, first ensure you have Python 3 installed and available on your PATH. Then, run the run script like so:

```bash
# Run a single sample
./run.sh dxsubmit --input input/config/inputs.full.json --dxworkflow <DNAnexus workflow path> --dxoutput <DNAnexus output path> --dxpriority [low,normal,high]

# Run a batch of samples
./run.sh dxsubmit --input input/config/inputs.full.json --batch input/inputFiles/inputFiles.tsv --dxworkflow <DNAnexus workflow path> --dxoutput <DNAnexus output path> --dxpriority [low,normal,high]
```

### Input parameters

The following tables describe the main parameters for each of the pipeline modes.

##### PoN creation

| Parameter(s) | Type | Description | Optional/Required | Default | Example |
| ------------ | ---- | ----------- | ----------------- | ------- | ------- |
| intervals | File | A Picard-formatted intervals list file | Optional, but recommended to reduce computation time and cost. | N/A | https://storage.googleapis.com/kccg-somvar-data/exome_evaluation_regions.v1.interval_list |
| ref_fasta, ref_fai, ref_dict | File | Reference genome FASTA, index, and dictionary | Required | N/A | https://storage.googleapis.com/kccg-somvar-data/Homo_sapiens_assembly38.fasta, https://storage.googleapis.com/kccg-somvar-data/Homo_sapiens_assembly38.fasta.fai, https://storage.googleapis.com/kccg-somvar-data/Homo_sapiens_assembly38.dict |
| bam_bai_list | File | Input TSV file containing the locations of the input BAMs and their index files; the file should have two columns, the first listing the BAM files and the second listing the indexes; one line per sample | Optional; if not supplied, `vcf_idx_list` must be supplied; both `bam_bai_list` and `vcf_idx_list` can be specified at the same time | N/A |  |
| vcf_idx_list | File | Input TSV file containing the locations of the pre-called input VCFs and their index files; the file should have two columns, the first listing the VCF files and the second listing the indexes; one line per sample | Optional; if not supplied, `bam_bai_list` must be supplied; both `bam_bai_list` and `vcf_idx_list` can be specified at the same time | N/A |  |
| scatter_count | Integer | The number of parallel shards to split the Mutect2 somatic variant calling job into | Optional | 10 |  |
| gnomad, gnomad_idx | File | Germline reference VCF and index containing common and rare variant population allele frequencies | Optional, but recommended | N/A | https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz, https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi |
| m2_extra_args | String | Extra arguments to pass to Mutect2 | Optional | N/A | "--pcr-indel-model NONE --downsampling-stride 20 --max-reads-per-alignment-start 6 --max-suspicious-reads-per-alignment-start 6" |
| compress | Boolean | Specify whether to compress output VCF files | Optional | false | true/false |
| samtools_docker | String | URI for a Docker image for running samtools commands | Optional | "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/vep-loftee@sha256:c95b78bacef4c8d3770642138e6f28998a5034cfad3fbef5451d2303c8c795d3" |  |
| gatk_docker | String | URI for a Docker image for running GATK commands | Optional | "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/gatk@sha256:0359ae4f32f2f541ca86a8cd30ef730bbaf8c306b9d53d2d520262d3e84b3b2b" |  |
| gatk_override | File | Path to an optional GATK .jar file to override the default container GATK version | Optional | N/A |  |
| preemptible | Integer | Number of times to allow a Google Cloud job to be pre-empted before running on a non-pre-embitble machine | Optional | 2 |  |
| max_retires | Integer | Number of times to allow a Google Cloud jobs to fail (not due to pre-empting) before giving up | Optional | 2 |  |
| small_task_cpu, small_task_mem, small_task_disk | Integer | Number of CPUs, amount of memory (MB), and amount of disk space (GB) to give most tasks | Optional | 4, 4000, 100 |  |
| command_mem_padding | Integer | Amount of memory (MB) to reserve for the Java runtime; the amount of memory given to run the command will be the total task memory minus the command_mem_padding value | Optional | 1000 |  |
| boot_disk_size | Integer | Amount of boot disk space to request when starting a job on Google Cloud | Optional | 12 |  |
| c2b_mem | Integer | Amount of memory (MB) to give to CramToBam | Optional | N/A |  |
| m2_mem, m2_cpu | Integer | Amount of memory (MB) and number of CPUs to use when running Mutect2 | Optional | 5000, 4 |  |
| emergency_extra_disk | Integer | Extra disk space to give to jobs in case jobs are failing due to running out of disk space | Optional | N/A |  |
| small_input_to_output_multiplier | Float | Disk space multiplier for small jobs; used to account for output size being larger than input size | Optional | 2.0 |  |
| cram_to_bam_multiplier | Float | Disk space multiplier for CramToBam; used to account for output size being larger than input size | Optional | N/A |  |
| m2_only | Boolean | Only run Mutect2, don't create the PoN | Optional | false | true/false |
| create_pon_extra_args | String | Additional arguments to pass to CreateSomaticPanelOfNormals | Optional | N/A |  |
| pon_name | String | Output name of the PoN VCF | Required | N/A |  |
| min_contig_size | Integer | Specify the minimum contig size when splitting genomic intervals | Optional | 1000000 |  |
| create_panel_scatter_count | Integer | The number of parallel shards to split the CreateSomaticPanelOfNormals somatic variant calling job into | Optional | 24 |  |

##### Full pipeline

| Parameter(s) | Type | Description | Optional/Required | Default | Example |
| ------------ | ---- | ----------- | ----------------- | ------- | ------- |
| intervals | File | A Picard-formatted intervals list file | Optional, but recommended to reduce computation time and cost. | N/A | https://storage.googleapis.com/kccg-somvar-data/exome_evaluation_regions.v1.interval_list |
| ref_fasta, ref_fai, ref_dict | File | Reference genome FASTA, index, and dictionary | Required | N/A | https://storage.googleapis.com/kccg-somvar-data/Homo_sapiens_assembly38.fasta, https://storage.googleapis.com/kccg-somvar-data/Homo_sapiens_assembly38.fasta.fai, https://storage.googleapis.com/kccg-somvar-data/Homo_sapiens_assembly38.dict |
| tumor_reads, tumor_reads_index | File | **(Single-sample mode only)** BAM alignments and BAI index for tumor sample | Required for single sample mode | N/A |  |
| normal_reads, normal_reads_index | File | **(Single-sample mode only)** BAM alignments and BAI index for matched normal sample | Optional; single sample mode only | N/A |  |
| pair_list | File | **(Batch mode only)** TSV file listing all tumor (and optional matched normal) sample BAM and BAI files | Required for batch mode | N/A |  |
| pon, pon_idx | File | Panel of normals VCF and VCF index files | Optional, but recommended | N/A | https://storage.googleapis.com/kccg-somvar-data/1000g_pon.hg38.vcf.gz, https://storage.googleapis.com/kccg-somvar-data/1000g_pon.hg38.vcf.gz.tbi |
| scatter_count | Integer | The number of parallel shards to split the Mutect2 somatic variant calling job into | Optional | 10 |  |
| gnomad, gnomad_idx | File | Germline reference VCF and index containing common and rare variant population allele frequencies | Optional, but recommended | N/A | https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz, https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi |
| variants_for_contamination, variants_for_contamination_idx | File | VCF and index containing common variants and allele frequencies for calculating contamination | Optional | N/A | https://storage.googleapis.com/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz, https://storage.googleapis.com/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz.tbi |
| realignment_index_bundle | File | BWA-mem index image to pass to FilterAlignmentArtifacts; causes FilterAlignmentArtifacts to run (experimental and not recommended) | Optional (not recommended) | N/A |  |
| realignment_extra_args | String | Extra arguments to pass to FilterAlignmentArtifacts | Optional | N/A |  |
| run_orientation_bias_mixture_model_filter | Boolean | Specify whether to run LearnReadOrientationModel | Optional | false | true/false |
| m2_extra_args | String | Extra arguments to pass to Mutect2 | Optional | N/A | "--pcr-indel-model NONE --downsampling-stride 20 --max-reads-per-alignment-start 6 --max-suspicious-reads-per-alignment-start 6" |
| m2_extra_filtering_args | String | Extra arguments to pass to FilterMutectCalls | Optional | N/A |  |
| split_intervals_extra_args | String | Extra arguments to pass to SplitIntervals | Optional | N/A |  |
| make_bamout | Boolean | Specify whether to gather and merge the Mutect2 output BAM files | Optional | false | true/false |
| compress_vcfs | Boolean | Specify whether to compress output VCF files | Optional | false | true/false |
| gga_vcf, gga_vcf_idx | File | Specify a set of alleles to force-call, regardless of evidence | Optional | N/A |  |
| vep | Boolean | Specifies whether or not to run VEP annotation | Optional | true | true/false |
| loftee | Boolean | Specifies whether or not to use LOFTEE annotations when running VEP; requires `vep` to be `true` | Optional | true | true/false |
| vep_docker | String | URI for a Docker image for running VEP | Optional | "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/vep@sha256:bc6a74bf271adb1484ea769660c7b69f5eea033d3ba2e2947988e6c5f034f221" |  |
| loftee_docker | String | URI for a Docker image for running LOFTEE | Optional | "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/vep-loftee@sha256:c95b78bacef4c8d3770642138e6f28998a5034cfad3fbef5451d2303c8c795d3" |  |
| vep_species | String | Species name to supply to the VEP `--species` parameter | Optional | "homo_sapiens" |  |
| vep_assembly | String | Assembly name to supply to the VEP `--assembly` parameter | Optional | "GRCh38" |  |
| vep_cache_archive | File | TAR.GZ archive file containing a VEP cache for annotating variants offline | Required for running vep (`vep` is set to `true`) | N/A | http://ftp.ensembl.org/pub/release-103/variation/vep/homo_sapiens_vep_103_GRCh38.tar.gz |
| vep_loftee_ancestor_fa, vep_loftee_ancestor_fai, vep_loftee_ancestor_gzi | File | FASTA file and index for running VEP + LOFTEE | Required for running VEP + Loftee (`vep` and `loftee` are set to `true`) | N/A | https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz, https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz.fai, https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz.gzi |
| vep_loftee_conservation_sql | File | PhyloCSF database for conservation filters in VEP + LOFTEE | Required for running VEP + Loftee (`vep` and `loftee` are set to `true`) | N/A | https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/loftee.sql.gz |
| annovar | Boolean | Specifies whether or not to run annovar annotation | Optional | false | true/false |
| annovar_docker | String | URI for a Docker image for running annovar | Optional | "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/perl@sha256:1f35086e2ff48dace3b3edeaa2ad1faf1e44c0612e00f00ea0fc1830b576a261" |  |
| annovar_assembly | String | Assembly name to supply to the annovar `-buildver` parameter | Optional | "hg38" |  |
| annovar_protocols, annovar_operations | List of annovar protocols and their respective operation codes to supply to the annovar parameters `-protocol` and `-operation`, respectively; must be comma-delimited and of the same length | Optional | `annovar_protocols`: "cosmic70"; `annovar_operations`: "f" |  |  |
| annovar_archive | File | TAR.GZ archive file containing the necessary files to run ANNOVAR | Required for running Annovar or CHIP detection (`annovar` is set to `true` OR `run_chip_detection` is set to `true`) | N/A | https://storage.cloud.google.com/kccg-somvar-data/annovar_files.tar.gz |
| run_chip_detection | Boolean | Specifies whether or not to run the CHIP detection script | Optional | true | true/false |
| treat_missing_as_rare | Boolean | (CHIP detection stage) When gathering allele frequencies from gnomAD, interpret missing variants as having an AF=0 | Optional | true | true/false |
| whitelist_genome | Boolean | (CHIP detection stage) Use gnomAD genome annotations as well as exome annotations | Optional | true | true/false |
| whitelist_use_ensembl_annotation | Boolean | (CHIP detection stage) Use Ensembl variant annotations rather than the default RefSeq annotations | Optional | false | true/false |
| run_chip_on_unannotated_vcf | Boolean | (CHIP detection stage) Don't use the annovar/VEP annotated VCF as input to CHIP; this won't change the CHIP detection algorithm, but will reduce the input file size to CHIP; annotated VCFs will still be generated | Optional | false | true/false |
| gnomad_pop | String | The gnomAD population code to use for gathering allele frequencies | Optional | "AF" (total AF, not sub-population-specific) |  |
| whitelist_docker | String | URI for a Docker image for running the CHIP detection script | Optional | "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/whitelist_filter@sha256:9cd77186c23a0b256a0928c5a4087b8378c234cb0754f35557cf9ec6d4aa544d" |  |
| whitelist_filter_archive | File | TAR.GZ archive file containing the necessary files to run the CHIP detection and variant whitelisting | Required for running CHIP detection (`run_chip_detection` is set to `true`) | N/A | https://storage.cloud.google.com/kccg-somvar-data/whitelist_filter_files.tar.gz |
| samtools_docker | String | URI for a Docker image for running samtools commands | Optional | "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/vep-loftee@sha256:c95b78bacef4c8d3770642138e6f28998a5034cfad3fbef5451d2303c8c795d3" |  |
| gatk_docker | String | URI for a Docker image for running GATK commands | Optional | "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/gatk@sha256:0359ae4f32f2f541ca86a8cd30ef730bbaf8c306b9d53d2d520262d3e84b3b2b" |  |
| gatk_override | File | Path to an optional GATK .jar file to override the default container GATK version | Optional | N/A |  |
| preemptible | Integer | Number of times to allow a Google Cloud job to be pre-empted before running on a non-pre-embitble machine | Optional | 2 |  |
| max_retires | Integer | Number of times to allow a Google Cloud jobs to fail (not due to pre-empting) before giving up | Optional | 2 |  |
| small_task_cpu, small_task_mem, small_task_disk | Integer | Number of CPUs, amount of memory (MB), and amount of disk space (GB) to give most tasks | Optional | 4, 4000, 100 |  |
| command_mem_padding | Integer | Amount of memory (MB) to reserve for the Java runtime; the amount of memory given to run the command will be the total task memory minus the command_mem_padding value | Optional | 1000 |  |
| boot_disk_size | Integer | Amount of boot disk space to request when starting a job on Google Cloud | Optional | 12 |  |
| c2b_mem | Integer | Amount of memory (MB) to give to CramToBam | Optional | N/A |  |
| m2_mem, m2_cpu | Integer | Amount of memory (MB) and number of CPUs to use when running Mutect2 | Optional | 5000, 4 |  |
| learn_read_orientation_mem | Integer | Amount of memory (MB) to give to LearnReadOrientationModel | Optional | 5000 |  |
| filter_alignment_artifacts_mem | Integer | Amount of memory (MB) to give to FilterAlignmentArtifacts | Optional | 5000 |  |
| vep_mem, vep_cpu, vep_tmp_disk | Integer | Amount of memory (MB), number of CPUs, and amount of temporary disk space (GB) to give to VEP | Optional | 32000, 1, 100 |  |
| annovar_mem_mb, annovar_disk, annovar_tmp_disk | Integer | Amount of memory (MB), disk space (GB), and temporary disk space (GB) to give to annovar | Optional | 4000, 100, 200 |  |
| whitelist_mem_mb, whitelist_disk | Integer | Amount of memory (MB) and disk space (GB) to give to the CHIP detection stage | Optional | 10000, 300 |  |
| emergency_extra_disk | Integer | Extra disk space to give to jobs in case jobs are failing due to running out of disk space | Optional | 0 |  |
| small_input_to_output_multiplier | Float | Disk space multiplier for small jobs; used to account for output size being larger than input size | Optional | 2.0 |  |
| large_input_to_output_multiplier | Float | Disk space multiplier for large jobs (currently only MergeBamOuts); used to account for output size being larger than input size | Optional | 2.25 |  |
| cram_to_bam_multiplier | Float | Disk space multiplier for CramToBam | Optional | N/A |  |

##### Annovar only

| Parameter(s) | Type | Description | Optional/Required | Default | Example |
| ------------ | ---- | ----------- | ----------------- | ------- | ------- |
| input_vcf | File | **(Single-sample mode only)** VCF file for annotation | Required for single sample mode | N/A |  |
| input_vcf_list | File | **(Batch mode only)** TSV file containing list of VCFs to annotate (one per line) | Required for batch mode | N/A |  |
| annovar_docker | String | URI for a Docker image for running annovar | Optional | "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/perl@sha256:1f35086e2ff48dace3b3edeaa2ad1faf1e44c0612e00f00ea0fc1830b576a261" |  |
| ref_name | String | Assembly name to supply to the annovar `-buildver` parameter | Optional | "hg38" |  |
| annovar_protocols, annovar_operations | List of annovar protocols and their respective operation codes to supply to the annovar parameters `-protocol` and `-operation`, respectively; must be comma-delimited and of the same length | Optional | `annovar_protocols`: "cosmic70"; `annovar_operations`: "f" |  |  |
| annovar_archive | File | TAR.GZ archive file containing the necessary files to run ANNOVAR | Required for running Annovar or CHIP detection (`annovar` is set to `true` OR `run_chip_detection` is set to `true`) | N/A | https://storage.cloud.google.com/kccg-somvar-data/annovar_files.tar.gz |
| preemptible | Integer | Number of times to allow a Google Cloud job to be pre-empted before running on a non-pre-embitble machine | Optional | 2 |  |
| max_retires | Integer | Number of times to allow a Google Cloud jobs to fail (not due to pre-empting) before giving up | Optional | 2 |  |
| small_task_cpu, small_task_mem, small_task_disk | Integer | Number of CPUs, amount of memory (MB), and amount of disk space (GB) to give most tasks | Optional | 4, 4000, 100 |  |
| command_mem_padding | Integer | Amount of memory (MB) to reserve for the Java runtime; the amount of memory given to run the command will be the total task memory minus the command_mem_padding value | Optional | 1000 |  |
| boot_disk_size | Integer | Amount of boot disk space to request when starting a job on Google Cloud | Optional | 12 |  |
| annovar_cpu, annovar_mem_mb, annovar_disk, annovar_tmp_disk | Integer | Number of CPUs, amount of memory (MB), disk space (GB), and temporary disk space (GB) to give to annovar | Optional | 1, 4000, 100, 200 |  |
| emergency_extra_disk | Integer | Extra disk space to give to jobs in case jobs are failing due to running out of disk space | Optional | 0 |  |

##### VEP only

| Parameter(s) | Type | Description | Optional/Required | Default | Example |
| ------------ | ---- | ----------- | ----------------- | ------- | ------- |
| input_vcf, input_vcf_idx | File | **(Single-sample mode only)** VCF file and its index for annotation | Required for single sample mode | N/A |  |
| input_vcf_list | File | **(Batch mode only)** TSV file containing list of VCFs and their index files to annotate (one per line; first column for VCF files, second column for index files) | Required for batch mode | N/A |  |
| ref_fasta | File | Reference genome FASTA | Required | N/A | https://storage.googleapis.com/kccg-somvar-data/Homo_sapiens_assembly38.fasta, https://storage.googleapis.com/kccg-somvar-data/Homo_sapiens_assembly38.fasta.fai, https://storage.googleapis.com/kccg-somvar-data/Homo_sapiens_assembly38.dict |
| loftee | Boolean | Specifies whether or not to use LOFTEE annotations | Optional | true | true/false |
| vep_docker | String | URI for a Docker image for running VEP | Optional | "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/vep@sha256:bc6a74bf271adb1484ea769660c7b69f5eea033d3ba2e2947988e6c5f034f221" |  |
| loftee_docker | String | URI for a Docker image for running LOFTEE | Optional | "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/vep-loftee@sha256:c95b78bacef4c8d3770642138e6f28998a5034cfad3fbef5451d2303c8c795d3" |  |
| vep_species | String | Species name to supply to the VEP `--species` parameter | Optional | "homo_sapiens" |  |
| vep_assembly | String | Assembly name to supply to the VEP `--assembly` parameter | Optional | "GRCh38" |  |
| vep_cache_archive | File | TAR.GZ archive file containing a VEP cache for annotating variants offline | Required for running vep (`vep` is set to `true`) | N/A | http://ftp.ensembl.org/pub/release-103/variation/vep/homo_sapiens_vep_103_GRCh38.tar.gz |
| vep_loftee_ancestor_fa, vep_loftee_ancestor_fai, vep_loftee_ancestor_gzi | File | FASTA file and index for running VEP + LOFTEE | Required for running VEP + Loftee (`vep` and `loftee` are set to `true`) | N/A | https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz, https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz.fai, https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz.gzi |
| vep_loftee_conservation_sql | File | PhyloCSF database for conservation filters in VEP + LOFTEE | Required for running VEP + Loftee (`vep` and `loftee` are set to `true`) | N/A | https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/loftee.sql.gz |
| preemptible | Integer | Number of times to allow a Google Cloud job to be pre-empted before running on a non-pre-embitble machine | Optional | 2 |  |
| max_retires | Integer | Number of times to allow a Google Cloud jobs to fail (not due to pre-empting) before giving up | Optional | 2 |  |
| small_task_cpu, small_task_mem, small_task_disk | Integer | Number of CPUs, amount of memory (MB), and amount of disk space (GB) to give most tasks | Optional | 4, 4000, 100 |  |
| command_mem_padding | Integer | Amount of memory (MB) to reserve for the Java runtime; the amount of memory given to run the command will be the total task memory minus the command_mem_padding value | Optional | 1000 |  |
| boot_disk_size | Integer | Amount of boot disk space to request when starting a job on Google Cloud | Optional | 12 |  |
| vep_mem, vep_cpu, vep_tmp_disk | Integer | Amount of memory (MB), number of CPUs, and amount of temporary disk space (GB) to give to VEP | Optional | 32000, 1, 100 |  |
| emergency_extra_disk | Integer | Extra disk space to give to jobs in case jobs are failing due to running out of disk space | Optional | 0 |  |

##### CHIP only

| Parameter(s) | Type | Description | Optional/Required | Default | Example |
| ------------ | ---- | ----------- | ----------------- | ------- | ------- |
| input_vcf | File | **(Single-sample mode only)** VCF file for input into CHIP detection | Required for single sample mode | N/A |  |
| tumor_sample_name | String | **(Single-sample mode only)** Name of tumor sample as recorded in the VCF header | Required for single sample mode | N/A |  |
| input_vcf_list | File | **(Batch mode only)** TSV file containing list of VCFs and the matched tumor sample name (as recorded in the VCF header) for input into CHIP detection (one per line; first column for sample names, second column for VCF files) | Required for batch mode | N/A |  |
| ref_fasta | File | Reference genome FASTA | Required | N/A | https://storage.googleapis.com/kccg-somvar-data/Homo_sapiens_assembly38.fasta, https://storage.googleapis.com/kccg-somvar-data/Homo_sapiens_assembly38.fasta.fai, https://storage.googleapis.com/kccg-somvar-data/Homo_sapiens_assembly38.dict |
| annovar_docker | String | URI for a Docker image for running annovar (part of the CHIP detection pipeline) | Optional | "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/perl@sha256:1f35086e2ff48dace3b3edeaa2ad1faf1e44c0612e00f00ea0fc1830b576a261" |  |
| ref_name | String | Reference assembly name to supply to the annovar `-buildver` parameter and to the CHIP detection script | Optional | "hg38" |  |
| annovar_archive | File | TAR.GZ archive file containing the necessary files to run ANNOVAR | Required for running Annovar or CHIP detection (`annovar` is set to `true` OR `run_chip_detection` is set to `true`) | N/A | https://storage.cloud.google.com/kccg-somvar-data/annovar_files.tar.gz |
| treat_missing_as_rare | Boolean | (CHIP detection stage) When gathering allele frequencies from gnomAD, interpret missing variants as having an AF=0 | Optional | true | true/false |
| whitelist_genome | Boolean | (CHIP detection stage) Use gnomAD genome annotations as well as exome annotations | Optional | true | true/false |
| whitelist_use_ensembl_annotation | Boolean | (CHIP detection stage) Use Ensembl variant annotations rather than the default RefSeq annotations | Optional | false | true/false |
| gnomad_pop | String | The gnomAD population code to use for gathering allele frequencies | Optional | "AF" (total AF, not sub-population-specific) |  |
| whitelist_docker | String | URI for a Docker image for running the CHIP detection script | Optional | "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/whitelist_filter@sha256:9cd77186c23a0b256a0928c5a4087b8378c234cb0754f35557cf9ec6d4aa544d" |  |
| whitelist_archive, whitelist_archive | File | TAR.GZ archive file containing the necessary files to run the CHIP detection and variant whitelisting | Required for running CHIP detection (`run_chip_detection` is set to `true`) | N/A | https://storage.cloud.google.com/kccg-somvar-data/whitelist_filter_files.tar.gz |
| preemptible | Integer | Number of times to allow a Google Cloud job to be pre-empted before running on a non-pre-embitble machine | Optional | 2 |  |
| max_retires | Integer | Number of times to allow a Google Cloud jobs to fail (not due to pre-empting) before giving up | Optional | 2 |  |
| small_task_cpu, small_task_mem, small_task_disk | Integer | Number of CPUs, amount of memory (MB), and amount of disk space (GB) to give most tasks | Optional | 4, 4000, 100 |  |
| command_mem_padding | Integer | Amount of memory (MB) to reserve for the Java runtime; the amount of memory given to run the command will be the total task memory minus the command_mem_padding value | Optional | 1000 |  |
| boot_disk_size | Integer | Amount of boot disk space to request when starting a job on Google Cloud | Optional | 12 |  |
| annovar_cpu, annovar_mem_mb, annovar_disk, annovar_tmp_disk | Integer | Number of CPUs, amount of memory (MB), disk space (GB), and temporary disk space (GB) to give to annovar | Optional | 1, 4000, 100, 200 |  |
| whitelist_cpu, whitelist_mem_mb, whitelist_disk | Integer | Number of CPUs, amount of memory (MB) and disk space (GB) to give to the CHIP detection stage | Optional | 1, 10000, 300 |  |
| emergency_extra_disk | Integer | Extra disk space to give to jobs in case jobs are failing due to running out of disk space | Optional | 0 |  |

