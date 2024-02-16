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

struct AnnotationRuntime {
    Int max_retries
    Int preemptible
    Int boot_disk_size
}

workflow CHIPCohort {
    input {
        # Inputs
        File vcf_idx_list
        Float prevalence_threshold = 0.01
        
        # CHIP Annotation settings
        String chip_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/chip_annotation:latest"

        # Runtime options
        String bcftools_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/chip_pre_post_filter:latest"
        Int? preemptible
        Int? max_retries
        Int small_task_cpu = 1
        Int small_task_mem = 4000
        Int small_task_disk = 50
        Int boot_disk_size = 12
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
        call StripFilterInfo {
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

    call MergeStrippedVcfsAndSplitMultiAllelics {
        input:
            input_vcf_idx_pairs = StripFilterInfo.stripped_vcf_idx_pairs,
            cpu = small_task_cpu,
            mem_mb = small_task_mem,
            disk = small_task_disk,
            docker = bcftools_docker,
            runtime_params = standard_runtime
    }

    call CHIPCohortFilter {
        input:
            input_vcf = MergeStrippedVcfsAndSplitMultiAllelics.stripped_split_vcf,
            input_vcf_idx = MergeStrippedVcfsAndSplitMultiAllelics.stripped_split_vcf_idx,
            prevalence_threshold = prevalence_threshold,
            cpu = small_task_cpu,
            mem_mb = small_task_mem,
            disk = small_task_disk,
            docker = chip_docker,
            runtime_params = standard_runtime
    }

    scatter(row in vcf_idx_pairs) {
        File in_vcf = row[0]
        File in_vcf_idx = row[1]
        call AnnotateVcf {
            input:
                input_vcf = in_vcf,
                input_vcf_idx = in_vcf_idx,
                annotation_vcf = CHIPCohortFilter.annotation_vcf,
                annotation_vcf_idx = CHIPCohortFilter.annotation_vcf_idx,
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

task StripFilterInfo {
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
        bcftools annotate -x FMT,INFO,FILTER ~{input_vcf} | bgzip -c > ~{vcf_basename}.stripped.vcf.gz
        tabix -s1 -b2 -e2 ~{vcf_basename}.stripped.vcf.gz
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
        Array[File] stripped_vcf_idx_pairs = ["~{vcf_basename}.stripped.vcf.gz", "~{vcf_basename}.stripped.vcf.gz.tbi"]
    }
}

task MergeStrippedVcfsAndSplitMultiAllelics {
    input {
        Array[Array[File]] input_vcf_idx_pairs
        Int cpu
        Int mem_mb
        Int disk
        String docker
        AnnotationRuntime runtime_params
    }

    Array[String] input_vcfs = transpose(input_vcf_idx_pairs)[0]

    command <<<
        # Use bcftools to merge the input VCF files
        bcftools merge -m all ~{sep=" " input_vcfs} | \
            bcftools norm -m -any -O z -o cohort.stripped.split.vcf.gz
        tabix -s1 -b2 -e2 cohort.stripped.split.vcf.gz
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
        File stripped_split_vcf = "cohort.stripped.split.vcf.gz"
        File stripped_split_vcf_idx = "cohort.stripped.split.vcf.gz.tbi"
    }
}

task CHIPCohortFilter {
    input {
        File input_vcf
        File input_vcf_idx
        Float prevalence_threshold = 0.01
        Int cpu
        Int mem_mb
        Int disk
        String docker
        AnnotationRuntime runtime_params
    }

    command <<<
        gunzip -c ~{input_vcf} > cohort.vcf
        annotate_chip_cohort \
            --input_vcf cohort.vcf \
            --prevalence_threshold ~{prevalence_threshold}

        bcftools sort cohort.annotated.vcf > cohort.annotated.sorted.vcf
        bgzip -c cohort.annotated.sorted.vcf > cohort.annotated.sorted.vcf.gz
        tabix -s1 -b2 -e2 cohort.annotated.sorted.vcf.gz

        bcftools norm -m +any -O z -o cohort.annotated.merged.sorted.vcf.gz cohort.annotated.sorted.vcf.gz
        tabix -s1 -b2 -e2 cohort.annotated.merged.sorted.vcf.gz
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
        File annotation_vcf = "cohort.annotated.merged.sorted.vcf.gz"
        File annotation_vcf_idx = "cohort.annotated.merged.sorted.vcf.gz.tbi"
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
        bcftools annotate \
            -a ~{annotation_vcf} \
            -c "=FILTER,+INFO" \
            -o ~{vcf_basename}.chip.cohort_filter.vcf.gz \
            ~{input_vcf}
        tabix -s1 -b2 -e2 ~{vcf_basename}.chip.cohort_filter.vcf.gz
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
        Array[File] sites_only_vcf_idx_pairs = ["~{vcf_basename}.chip.cohort_filter.vcf.gz", "~{vcf_basename}.chip.cohort_filter.vcf.gz.tbi"]
    }
}