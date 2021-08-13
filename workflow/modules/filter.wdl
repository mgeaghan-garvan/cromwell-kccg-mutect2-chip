version 1.0

import "runtime.wdl" as RT

task Filter {
    input {
      File? intervals
      File ref_fasta
      File ref_fai
      File ref_dict
      File unfiltered_vcf
      File unfiltered_vcf_idx
      String output_name
      Boolean compress
      File? mutect_stats
      File? artifact_priors_tar_gz
      File? contamination_table
      File? maf_segments
      String? m2_extra_filtering_args

      Runtime runtime_params
      Int? disk_space
    }

    String output_vcf = output_name + if compress then ".vcf.gz" else ".vcf"
    String output_vcf_idx = output_vcf + if compress then ".tbi" else ".idx"

    parameter_meta{
      ref_fasta: {localization_optional: true}
      ref_fai: {localization_optional: true}
      ref_dict: {localization_optional: true}
    }

    command {
        set -e

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}

        gatk --java-options "-Xmx~{runtime_params.command_mem}m" FilterMutectCalls -V ~{unfiltered_vcf} \
            -R ~{ref_fasta} \
            -O ~{output_vcf} \
            ~{"--contamination-table " + contamination_table} \
            ~{"--tumor-segmentation " + maf_segments} \
            ~{"--ob-priors " + artifact_priors_tar_gz} \
            ~{"-stats " + mutect_stats} \
            --filtering-stats filtering.stats \
            ~{m2_extra_filtering_args}
    }

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
        File filtered_vcf = "~{output_vcf}"
        File filtered_vcf_idx = "~{output_vcf_idx}"
        File filtering_stats = "filtering.stats"
    }
}

task FilterAlignmentArtifacts {
    input {
      File ref_fasta
      File ref_fai
      File ref_dict
      File input_vcf
      File input_vcf_idx
      File bam
      File bai
      String output_name
      Boolean compress
      File realignment_index_bundle
      String? realignment_extra_args
      Runtime runtime_params
      Int mem
    }

    String output_vcf = output_name + if compress then ".vcf.gz" else ".vcf"
    String output_vcf_idx = output_vcf +  if compress then ".tbi" else ".idx"

    Int machine_mem = mem
    Int command_mem = machine_mem - 500

    parameter_meta{
      ref_fasta: {localization_optional: true}
      ref_fai: {localization_optional: true}
      ref_dict: {localization_optional: true}
      input_vcf: {localization_optional: true}
      input_vcf_idx: {localization_optional: true}
      bam: {localization_optional: true}
      bai: {localization_optional: true}
    }

    command {
        set -e

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}

        gatk --java-options "-Xmx~{command_mem}m" FilterAlignmentArtifacts \
            -R ~{ref_fasta} \
            -V ~{input_vcf} \
            -I ~{bam} \
            --bwa-mem-index-image ~{realignment_index_bundle} \
            ~{realignment_extra_args} \
            -O ~{output_vcf}
    }

    runtime {
        docker: runtime_params.gatk_docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: machine_mem + " MB"
        disks: "local-disk " + runtime_params.disk + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    output {
        File filtered_vcf = "~{output_vcf}"
        File filtered_vcf_idx = "~{output_vcf_idx}"
    }
}