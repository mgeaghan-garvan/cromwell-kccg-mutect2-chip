version 1.0

# =============================================================== #
# spliceai.wdl                                                    #
#                                                                 #
# This workflow adds SpliceAI annotations to a VCF.               #
#                                                                 #
# This pipeline has been developed for use by the Kinghorn        #
# Centre for Clinical Genomics and the Garvan Institute for       # 
# Medical Research.                                               #
#                                                                 #
# Author: Michael Geaghan (micgea)                                #
# Created: 2023/04/04                                             #
# =============================================================== #

struct SpliceAIRuntime {
    Int max_retries
    Int preemptible
    Int boot_disk_size
}

workflow SpliceAI {
    input {
        File input_vcf
        File ref_fasta
        File? spliceai_annotation_file
        String spliceai_annotation_string = 'grch38'
        Int spliceai_max_dist = 50
        Boolean spliceai_mask = false
        String spliceai_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/spliceai@sha256:617682496a3f475c69ccdfe593156b79dd1ba21e02481ed1d0d8b740f3422530"  # :v1.3.1
        Int spliceai_disk = 100
        Int spliceai_mem_mb = 16000
        Int spliceai_cpu = 4
        # runtime parameters
        Int? preemptible
        Int? max_retries
        Int boot_disk_size = 12
    }

    Int preemptible_or_default = select_first([preemptible, 2])
    Int max_retries_or_default = select_first([max_retries, 2])

    SpliceAIRuntime standard_runtime = {
        "max_retries": max_retries_or_default,
        "preemptible": preemptible_or_default,
        "boot_disk_size": boot_disk_size
    }

    call SpliceAI_task {
        input:
            input_vcf = input_vcf,
            ref_fasta = ref_fasta,
            annotation_file = spliceai_annotation_file,
            annotation_string = spliceai_annotation_string,
            max_dist = spliceai_max_dist,
            mask = spliceai_mask,
            spliceai_docker = spliceai_docker,
            disk_space = spliceai_disk,
            mem_mb = spliceai_mem_mb,
            cpu = spliceai_cpu,
            runtime_params = standard_runtime
    }

    output {
        File spliceai_output_vcf = SpliceAI_task.spliceai_output_vcf
    }
}

task SpliceAI_task {
    input {
      File input_vcf
      File ref_fasta
      File? annotation_file
      String annotation_string = 'grch38'
      Int max_dist = 50
      Boolean mask = false
      String spliceai_docker
      Int disk_space
      Int mem_mb
      Int cpu
      SpliceAIRuntime runtime_params
    }

    String input_basename = basename(basename(input_vcf, ".gz"), ".vcf")
    String annotation_param = select_first([annotation_file, annotation_string])
    Int mask_val = if mask then 1 else 0

    command {
      spliceai \
        -I ~{input_vcf} \
        -O ~{input_basename}.spliceai.vcf \
        -R ~{ref_fasta} \
        -A ~{annotation_param} \
        -D ~{max_dist} \
        -M ~{mask_val}
    }

    runtime {
      docker: spliceai_docker
      bootDiskSizeGb: runtime_params.boot_disk_size
      memory: mem_mb + " MB"
      disks: "local-disk " + disk_space + " HDD"
      preemptible: runtime_params.preemptible
      maxRetries: runtime_params.max_retries
      cpu: cpu
    }

    output {
      File spliceai_output_vcf = "~{input_basename}.spliceai.vcf"
    }
}