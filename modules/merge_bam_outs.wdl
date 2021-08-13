version 1.0

import "./runtime.wdl" as RT

task MergeBamOuts {
    input {
      File ref_fasta
      File ref_fai
      File ref_dict
      Array[File]+ bam_outs
      String output_vcf_name
      Runtime runtime_params
      Int? disk_space   #override to request more disk than default small task params
    }

    command <<<
        # This command block assumes that there is at least one file in bam_outs.
        #  Do not call this task if len(bam_outs) == 0
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" GatherBamFiles \
            -I ~{sep=" -I " bam_outs} -O unsorted.out.bam -R ~{ref_fasta}

        # We must sort because adjacent scatters may have overlapping (padded) assembly regions, hence
        # overlapping bamouts

        gatk --java-options "-Xmx~{runtime_params.command_mem}m" SortSam -I unsorted.out.bam \
            -O ~{output_vcf_name}.out.bam \
            --SORT_ORDER coordinate -VALIDATION_STRINGENCY LENIENT
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" BuildBamIndex -I ~{output_vcf_name}.out.bam -VALIDATION_STRINGENCY LENIENT
    >>>

    runtime {
        docker: runtime_params.gatk_docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, runtime_params.disk]) + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    output {
        File merged_bam_out = "~{output_vcf_name}.out.bam"
        File merged_bam_out_index = "~{output_vcf_name}.out.bai"
    }
}