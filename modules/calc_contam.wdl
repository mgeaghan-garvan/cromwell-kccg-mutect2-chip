version 1.0

import "./runtime.wdl" as RT

task CalculateContamination {
    input {
      String? intervals
      File tumor_pileups
      File? normal_pileups
      Runtime runtime_params
    }

    command {
        set -e

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}

        gatk --java-options "-Xmx~{runtime_params.command_mem}m" CalculateContamination -I ~{tumor_pileups} \
        -O contamination.table --tumor-segmentation segments.table ~{"-matched " + normal_pileups}
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
        File contamination_table = "contamination.table"
        File maf_segments = "segments.table"
    }
}