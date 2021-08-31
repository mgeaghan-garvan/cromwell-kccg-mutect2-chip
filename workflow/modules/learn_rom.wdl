version 1.0

import "runtime.wdl" as RT

# Learning step of the orientation bias mixture model, which is the recommended orientation bias filter as of September 2018
task LearnReadOrientationModel {
    input {
      Array[File] f1r2_tar_gz
      Runtime runtime_params
      Int mem_mb = 5000
    }

    Int machine_mem = mem_mb
    Int command_mem = machine_mem - 500

    command {
        set -e
        export GATK_LOCAL_JAR=~{default="/gatk/gatk.jar" runtime_params.gatk_override}

        gatk --java-options "-Xmx~{command_mem}m" LearnReadOrientationModel \
            -I ~{sep=" -I " f1r2_tar_gz} \
            -O "artifact-priors.tar.gz"
    }

    runtime {
        docker: runtime_params.gatk_docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        mem_mb: machine_mem
        disks: "local-disk " + runtime_params.disk + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    output {
        File artifact_prior_table = "artifact-priors.tar.gz"
    }

}