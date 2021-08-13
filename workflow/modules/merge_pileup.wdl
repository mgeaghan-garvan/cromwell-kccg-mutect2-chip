version 1.0

import "runtime.wdl" as RT

task MergePileupSummaries {
    input {
      Array[File] input_tables
      String output_name
      File ref_dict
      Runtime runtime_params
    }

    command {
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}

        gatk --java-options "-Xmx~{runtime_params.command_mem}m" GatherPileupSummaries \
        --sequence-dictionary ~{ref_dict} \
        -I ~{sep=' -I ' input_tables} \
        -O ~{output_name}.tsv
    }

    runtime {
        docker: runtime_params.gatk_docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.machine_mem + " MB"
        disks: "local-disk " + runtime_params.disk + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    output {
        File merged_table = "~{output_name}.tsv"
    }
}