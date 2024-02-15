version 1.0

# =============================================================== #
# cromwell-kccg-mutect2-chip                                      #
#                                                                 #
# This workflow detects mutations associated with CHIP            #
# (clonal haematopoiesis of indeterminate potential) using the    #
# GATK4 Mutect2 somatic variant calling pipeline and the          #
# Cromwell workflow engine.                                       #
#                                                                 #
# This workflow has been adapted from the CHIP-detection-Mutect2  #
# workflow developed by Alex Bick. The Mutect2 somatic variant    #
# calling stage is adapted from the GATK best practices workflow, #
# and the CHIP detection stage is adapted from the                #
# Annovar Whitelist Filter developed by Charlie Condon.           #
#                                                                 #
# This pipeline has been developed for use by the Kinghorn        #
# Centre for Clinical Genomics and the Garvan Institute for       # 
# Medical Research.                                               #
#                                                                 #
# Author: Michael Geaghan (micgea)                                #
# Created: 2021/08/13                                             #
# =============================================================== #

import "vep.wdl" as VEP
import "annovar.wdl" as Annovar
import "spliceai.wdl" as SpliceAI
import "chip.wdl" as CHIP

workflow AnnotateCohort {
    input {
        # Inputs
        File vcf_idx_list
        File ref_fasta
        File ref_fai
        String cohort_name = "cohort"
        
        # VEP settings
        Boolean vep = false
        File? vep_cache_archive
        String vep_species = "homo_sapiens"
        String vep_assembly = "GRCh38"
        String vep_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/vep@sha256:bc6a74bf271adb1484ea769660c7b69f5eea033d3ba2e2947988e6c5f034f221"  # :release_103.1
        Boolean loftee = false
        File? vep_loftee_ancestor_fa
        File? vep_loftee_ancestor_fai
        File? vep_loftee_ancestor_gzi
        File? vep_loftee_conservation_sql
        String loftee_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/vep-loftee@sha256:c95b78bacef4c8d3770642138e6f28998a5034cfad3fbef5451d2303c8c795d3"  # :vep_103.1_loftee_1.0.3

        # Annovar settings
        Boolean annovar = true
        File? annovar_db_archive
        String annovar_assembly = "hg38"
        String annovar_protocols = "refGene"
        String annovar_operations = "g"
        String annovar_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/annovar@sha256:842e9f88dd39999ee2129aeb992e8eced10ac2a33642d4b34d0f0c0254aa5035"  # :5.34.0

        # SpliceAI settings
        Boolean spliceai = false
        File? spliceai_annotation_file
        String spliceai_annotation_string = 'grch38'
        Int spliceai_max_dist = 50
        Boolean spliceai_mask = false
        String spliceai_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/spliceai@sha256:617682496a3f475c69ccdfe593156b79dd1ba21e02481ed1d0d8b740f3422530"  # :v1.3.1

        # CHIP Annotation settings
        File chip_mutations_csv
        File somaticism_filter_transcripts
        Boolean run_chip_detection = true
        Boolean treat_missing_as_rare = true
        Boolean use_gnomad_genome = true
        Boolean use_ensembl_annotation = false
        String gnomad_pop = "AF"
        String chip_pre_post_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/chip_pre_post_filter:latest"
        String chip_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/chip_annotation:latest"

        # Runtime options
        String bcftools_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/vep-loftee@sha256:c95b78bacef4c8d3770642138e6f28998a5034cfad3fbef5451d2303c8c795d3"  # same as loftee_docker
        Int? preemptible
        Int? max_retries
        Int small_task_cpu = 4
        Int small_task_mem = 4000
        Int small_task_disk = 100
        Int command_mem_padding = 1000
        Int boot_disk_size = 12
        Int vep_mem = 32000
        Int vep_cpu = 1
        Int vep_tmp_disk = 100
        Int annovar_mem_mb = 4000
        Int annovar_disk = 100
        Int spliceai_disk = 100
        Int spliceai_mem_mb = 16000
        Int spliceai_cpu = 4
        Int chip_mem_mb = 10000
        Int chip_disk = 300

        # Use as a last resort to increase the disk given to every task in case of ill behaving data
        Int? emergency_extra_disk

        # These are multipliers to multipler inputs by to make sure we have enough disk to accommodate for possible output sizes
        # Large is for Bams/WGS vcfs
        # Small is for metrics/other vcfs
        Float large_input_to_output_multiplier = 2.25
        Float small_input_to_output_multiplier = 2.0
    }

    Int preemptible_or_default = select_first([preemptible, 2])
    Int max_retries_or_default = select_first([max_retries, 2])

    AnnotationRuntime standard_runtime = {
        "max_retries": max_retries_or_default,
        "preemptible": preemptible_or_default,
        "boot_disk_size": boot_disk_size
    }

    Array[Array[String]] vcf_idx_pairs = read_tsv(select_first([vcf_idx_list, ""]))

    scatter(row in vcf_idx_pairs) {
        File vcf = row[0]
        File vcf_idx = row[1]
        call VcfToSitesOnlyVcf {
            input:
                input_vcf = vcf,
                input_vcf_idx = vcf_idx,
                cpu = small_task_cpu,
                mem_mb = small_task_mem,
                disk = small_task_disk,
                docker = bcftools_docker,
                runtime_params = standard_runtime
        }
    }

    call MergeSitesOnlyVcfs {
        input:
            input_vcf_idx_pairs = VcfToSitesOnlyVcf.sites_only_vcf_idx_pairs,
            cohort_name = cohort_name,
            cpu = small_task_cpu,
            mem_mb = small_task_mem,
            disk = small_task_disk,
            docker = bcftools_docker,
            runtime_params = standard_runtime
    }

    # Optionally run VEP
    if (vep && defined(vep_cache_archive)) {
        File vep_cache_archive_file = select_first([vep_cache_archive, "CACHE_FILE_NOT_SUPPLIED"])
        call VEP.VEP as VEP_wf {
            input:
                input_vcf = MergeSitesOnlyVcfs.sites_only_vcf,
                input_vcf_idx = MergeSitesOnlyVcfs.sites_only_vcf_idx,
                vep_species = vep_species,
                vep_assembly = vep_assembly,
                vep_cache_archive = vep_cache_archive_file,
                ref_fasta = ref_fasta,
                vep_docker = vep_docker,
                loftee = loftee,
                vep_loftee_ancestor_fa = vep_loftee_ancestor_fa,
                vep_loftee_ancestor_fai = vep_loftee_ancestor_fai,
                vep_loftee_ancestor_gzi = vep_loftee_ancestor_gzi,
                vep_loftee_conservation_sql = vep_loftee_conservation_sql,
                loftee_docker = loftee_docker,
                vep_mem = vep_mem,
                vep_disk = small_task_disk,
                vep_tmp_disk = vep_tmp_disk,
                vep_cpu = vep_cpu,
                preemptible = preemptible,
                max_retries = max_retries,
                boot_disk_size = boot_disk_size,
                emergency_extra_disk = emergency_extra_disk
        }
    }

    File annovar_input_vcf = select_first([VEP_wf.out_vep_vcf, MergeSitesOnlyVcfs.sites_only_vcf])
    File annovar_db_archive_file = select_first([annovar_db_archive, "ANNOVAR_ARCHIVE_NOT_SUPPLIED"])

    # Optionally run Annovar
    if (annovar && defined(annovar_db_archive)) {
        call Annovar.Annovar as Annovar_wf {
            input:
                input_vcf = annovar_input_vcf,
                annovar_mem_mb = annovar_mem_mb,
                annovar_disk = annovar_disk,
                annovar_cpu = 1,
                annovar_docker = annovar_docker,
                annovar_db_archive = annovar_db_archive_file,
                ref_name = annovar_assembly,
                annovar_protocols = annovar_protocols,
                annovar_operations = annovar_operations,
                preemptible = preemptible,
                max_retries = max_retries,
                boot_disk_size = boot_disk_size
        }
    }

    File splieai_input_vcf = select_first([Annovar_wf.out_annovar_vcf, annovar_input_vcf])

    # Optionally run SpliceAI
    if (spliceai) {
        call SpliceAI.SpliceAI as SpliceAI_wf {
            input:
                input_vcf = splieai_input_vcf,
                ref_fasta = ref_fasta,
                spliceai_annotation_file = spliceai_annotation_file,
                spliceai_annotation_string = spliceai_annotation_string,
                spliceai_max_dist = spliceai_max_dist,
                spliceai_mask = spliceai_mask,
                spliceai_docker = spliceai_docker,
                spliceai_disk = spliceai_disk,
                spliceai_mem_mb = spliceai_mem_mb,
                spliceai_cpu = spliceai_cpu,
                preemptible = preemptible,
                max_retries = max_retries,
                boot_disk_size = boot_disk_size
        }
    }

    File annotation_vcf = select_first([SpliceAI_wf.spliceai_output_vcf, splieai_input_vcf])

    call IndexVcf {
        input:
            input_vcf = annotation_vcf,
            cpu = small_task_cpu,
            mem_mb = small_task_mem,
            disk = small_task_disk,
            docker = bcftools_docker,
            runtime_params = standard_runtime
    }

    File annotation_vcf_idx = IndexVcf.output_vcf_idx

    scatter(row in vcf_idx_pairs) {
        File in_vcf = row[0]
        File in_vcf_idx = row[1]
        call AnnotateVcf {
            input:
                input_vcf = in_vcf,
                input_vcf_idx = in_vcf_idx,
                annotation_vcf = annotation_vcf,
                annotation_vcf_idx = annotation_vcf_idx,
                cpu = small_task_cpu,
                mem_mb = small_task_mem,
                disk = small_task_disk,
                docker = bcftools_docker,
                runtime_params = standard_runtime
        }
    }

    output {
        Array[Array[File]] annotated_vcf_idx_pairs = AnnotateVcf.sites_only_vcf_idx_pairs
    }
}

task VcfToSitesOnlyVcf {
    input {
        File input_vcf
        File input_vcf_idx
        Int cpu
        Int mem_mb
        Int disk
        String docker
        AnnotationRuntime runtime_params
    }

    String vcf_basename = basename(basename(input_vcf, ".gz"), ".vcf")

    command <<<
        # Use bcftools to convert the input VCF to a sites-only VCF
        bcftools annotate -x FMT,INFO,FILTER ~{input_vcf} | cut -f 1-8 | bgzip -c > ~{vcf_basename}.sites_only.vcf.gz
        tabix -s1 -b2 -e2 ~{vcf_basename}.sites_only.vcf.gz
    >>>

    runtime {
        docker: docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: mem_mb + " MB"
        disks: "local-disk " + disk + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: cpu
    }

    output {
        Array[File] sites_only_vcf_idx_pairs = ["~{vcf_basename}.sites_only.vcf.gz", "~{vcf_basename}.sites_only.vcf.gz.tbi"]
    }
}

task MergeSitesOnlyVcfs {
    input {
        Array[Array[File]] input_vcf_idx_pairs
        String cohort_name
        Int cpu
        Int mem_mb
        Int disk
        String docker
        AnnotationRuntime runtime_params
    }

    Array[String] input_vcfs = transpose(input_vcf_idx_pairs)[0]

    command <<<
        # Use bcftools to merge the input VCF files
        bcftools concat -a -d exact -O z -o ~{cohort_name}.sites_only.vcf.gz ~{sep=" " input_vcfs}
        tabix -s1 -b2 -e2 ~{cohort_name}.sites_only.vcf.gz
    >>>

    runtime {
        docker: docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: mem_mb + " MB"
        disks: "local-disk " + disk + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: cpu
    }

    output {
        File sites_only_vcf = "~{cohort_name}.sites_only.vcf.gz"
        File sites_only_vcf_idx = "~{cohort_name}.sites_only.vcf.gz.tbi"
    }
}

task IndexVcf {
    input {
        File input_vcf
        Int cpu
        Int mem_mb
        Int disk
        String docker
        AnnotationRuntime runtime_params
    }

    command <<<
        tabix -s1 -b2 -e2 ~{input_vcf}
    >>>

    runtime {
        docker: docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: mem_mb + " MB"
        disks: "local-disk " + disk + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: cpu
    }

    output {
        File output_vcf_idx = input_vcf + ".tbi"
    }

}

task AnnotateVcf {
    input {
        File input_vcf
        File input_vcf_idx
        File annotation_vcf
        File annotation_vcf_idx
        Int cpu
        Int mem_mb
        Int disk
        String docker
        AnnotationRuntime runtime_params
    }

    String vcf_basename = basename(basename(input_vcf, ".gz"), ".vcf")

    command <<<
        # Use bcftools to annotate the input VCF with the annotation VCF
        bcftools annotate -a ~{annotation_vcf} -c CHROM,POS,REF,ALT,INFO ~{input_vcf} | bgzip -c > ~{vcf_basename}.annotated.vcf.gz
        tabix -s1 -b2 -e2 ~{vcf_basename}.annotated.vcf.gz
    >>>

    runtime {
        docker: docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: mem_mb + " MB"
        disks: "local-disk " + disk + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: cpu
    }

    output {
        Array[File] sites_only_vcf_idx_pairs = ["~{vcf_basename}.annotated.vcf.gz", "~{vcf_basename}.annotated.vcf.gz.tbi"]
    }
}