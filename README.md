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

##### Input parameters

The following tables describe the main parameters for each of the pipeline modes.

###### PoN creation

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

###### Full pipeline

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
| compress_vcfs, compress | Boolean | Specify whether to compress output VCF files | Optional | false | true/false |
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
| whitelist_docker | String | URI for a Docker image for running the CHIP detection script | Optional | "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/whitelist_filter@sha256:1f1f83f8241f40fbd1f21b19e2ccbdc184984fd9ec0b0a7bdfa97b8a73fed8a4" |  |
| whitelist_filter_archive, whitelist_archive | File | TAR.GZ archive file containing the necessary files to run the CHIP detection and variant whitelisting | Required for running CHIP detection (`run_chip_detection` is set to `true`) | N/A | https://storage.cloud.google.com/kccg-somvar-data/whitelist_filter_files.tar.gz |
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

###### Annovar only

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

###### VEP only

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

###### CHIP only

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
| whitelist_docker | String | URI for a Docker image for running the CHIP detection script | Optional | "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/whitelist_filter@sha256:1f1f83f8241f40fbd1f21b19e2ccbdc184984fd9ec0b0a7bdfa97b8a73fed8a4" |  |
| whitelist_archive, whitelist_archive | File | TAR.GZ archive file containing the necessary files to run the CHIP detection and variant whitelisting | Required for running CHIP detection (`run_chip_detection` is set to `true`) | N/A | https://storage.cloud.google.com/kccg-somvar-data/whitelist_filter_files.tar.gz |
| preemptible | Integer | Number of times to allow a Google Cloud job to be pre-empted before running on a non-pre-embitble machine | Optional | 2 |  |
| max_retires | Integer | Number of times to allow a Google Cloud jobs to fail (not due to pre-empting) before giving up | Optional | 2 |  |
| small_task_cpu, small_task_mem, small_task_disk | Integer | Number of CPUs, amount of memory (MB), and amount of disk space (GB) to give most tasks | Optional | 4, 4000, 100 |  |
| command_mem_padding | Integer | Amount of memory (MB) to reserve for the Java runtime; the amount of memory given to run the command will be the total task memory minus the command_mem_padding value | Optional | 1000 |  |
| boot_disk_size | Integer | Amount of boot disk space to request when starting a job on Google Cloud | Optional | 12 |  |
| annovar_cpu, annovar_mem_mb, annovar_disk, annovar_tmp_disk | Integer | Number of CPUs, amount of memory (MB), disk space (GB), and temporary disk space (GB) to give to annovar | Optional | 1, 4000, 100, 200 |  |
| whitelist_cpu, whitelist_mem_mb, whitelist_disk | Integer | Number of CPUs, amount of memory (MB) and disk space (GB) to give to the CHIP detection stage | Optional | 1, 10000, 300 |  |
| emergency_extra_disk | Integer | Extra disk space to give to jobs in case jobs are failing due to running out of disk space | Optional | 0 |  |

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
