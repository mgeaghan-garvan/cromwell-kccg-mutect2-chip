version 1.0

# ======================================================================================= #
# pon.wdl                                                                                 #
#                                                                                         #
# This workflow has been adapted from the mutect2_pon.wdl workflow developed by the       #
# Broad Institue. The base version of the original workflow can be found here:            #
# https://github.com/broadinstitute/gatk/blob/4.1.6.0/scripts/mutect2_wdl/mutect2_pon.wdl #
#                                                                                         #
# A few minor alterations have been made for integration into the larger CHIP workflow    #
#                                                                                         #
# This pipeline has been developed for use by the Kinghorn                                #
# Centre for Clinical Genomics and the Garvan Institute for                               # 
# Medical Research.                                                                       #
#                                                                                         #
# Author: Michael Geaghan (micgea)                                                        #
# Created: 2023/04/04                                                                     #
# ======================================================================================= #

#  Create a Mutect2 panel of normals
#
#  Description of inputs
#  intervals: genomic intervals
#  ref_fasta, ref_fai, ref_dict: reference genome, index, and dictionary
#  normal_bams: arrays of normal bams
#  scatter_count: number of parallel jobs when scattering over intervals
#  pon_name: the resulting panel of normals is {pon_name}.vcf
#  m2_extra_args: additional command line parameters for Mutect2.  This should not involve --max-mnp-distance,
#  which the wdl hard-codes to 0 because GenomicsDBImport can't handle MNPs

import "m2.wdl" as m2
import "call_u2af1.wdl" as U2AF1

workflow Mutect2_Panel {
    input {
        File? intervals
        File ref_fasta
        File ref_fai
        File ref_dict
        File? bam_bai_list
        File? vcf_idx_list
        Int scatter_count
        File gnomad
        File gnomad_idx
        String? m2_extra_args
        Boolean m2_only = false
        String? create_pon_extra_args
        Boolean? compress
        String pon_name

        Int? min_contig_size
        Int? create_panel_scatter_count

        # U2AF1 settings
        Boolean u2af1 = true
        File u2af1_regions_file
        String pileup_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/u2af1:latest"
        String merge_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/chip_pre_post_filter:latest"

        # Runtime options
        String gatk_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/gatk@sha256:0359ae4f32f2f541ca86a8cd30ef730bbaf8c306b9d53d2d520262d3e84b3b2b"  # :4.2.1.0
        File? gatk_override
        Int? preemptible
        Int? max_retries
        Int small_task_cpu = 4
        Int small_task_mem = 4000
        Int small_task_disk = 100
        Int command_mem_padding = 1000
        Int boot_disk_size = 12
        Int m2_mem = 5000
        Int m2_cpu = 4
        Int pon_mem = 5000

        # These are multipliers to multipler inputs by to make sure we have enough disk to accommodate for possible output sizes
        Float small_input_to_output_multiplier = 2.0
        Float cram_to_bam_multiplier = 6.0
    }

    Int contig_size = select_first([min_contig_size, 1000000])
    Int preemptible_or_default = select_first([preemptible, 2])
    Int max_retries_or_default = select_first([max_retries, 2])

    Int small_task_cpu_mult = if small_task_cpu > 1 then small_task_cpu - 1 else 1
    Int cmd_mem = small_task_mem - command_mem_padding

    Runtime standard_runtime = {
        "gatk_docker": gatk_docker,
        "gatk_override": gatk_override,
        "max_retries": max_retries_or_default,
        "preemptible": preemptible_or_default,
        "cpu": small_task_cpu,
        "machine_mem": small_task_mem,
        "command_mem": cmd_mem,
        "disk": small_task_disk,
        "boot_disk_size": boot_disk_size
    }

    Array[Array[String]] bam_bai_pairs = if (defined(bam_bai_list)) then read_tsv(select_first([bam_bai_list, ""])) else [[]]
    Array[Array[String]] vcf_idx_pairs = if (defined(vcf_idx_list)) then transpose(read_tsv(select_first([vcf_idx_list, ""]))) else [[], []]

    if (defined(bam_bai_list)) {
        # scatter (normal_bam in zip(normal_bams, normal_bais)) {
        scatter (row in bam_bai_pairs) {
            File normal_bam = row[0]
            File normal_bai = row[1]
            call m2.Mutect2 {
                input:
                    intervals = intervals,
                    ref_fasta = ref_fasta,
                    ref_fai = ref_fai,
                    ref_dict = ref_dict,
                    tumor_reads = normal_bam,
                    tumor_reads_index = normal_bai,
                    scatter_count = scatter_count,
                    m2_extra_args = select_first([m2_extra_args, ""]) + " --max-mnp-distance 0",
                    compress_vcfs = select_first([compress, false]),
                    gatk_docker = gatk_docker,
                    gatk_override = gatk_override,
                    preemptible = preemptible,
                    max_retries = max_retries,
                    small_task_cpu = small_task_cpu,
                    small_task_mem = small_task_mem,
                    small_task_disk = small_task_disk,
                    command_mem_padding = command_mem_padding,
                    boot_disk_size = boot_disk_size,
                    m2_mem = m2_mem,
                    m2_cpu = m2_cpu
            }

            if (u2af1 && defined(u2af1_regions_file)) {
                call U2AF1.CallU2AF1 as U2AF1_wf {
                    input:
                        tumor_reads = normal_bam,
                        tumor_reads_index = normal_bai,
                        mutect2_output_vcf = Mutect2.filtered_vcf,
                        mutect2_output_vcf_index = Mutect2.filtered_vcf_idx,
                        ref_fasta = ref_fasta,
                        ref_fai = ref_fai,
                        u2af1_regions_file = u2af1_regions_file,
                        pileup_docker = pileup_docker,
                        merge_docker = merge_docker,
                        preemptible = preemptible,
                        max_retries = max_retries,
                        cpu = small_task_cpu,
                        mem_mb = small_task_mem,
                        disk = small_task_disk,
                        boot_disk_size = boot_disk_size
                }
            }

            File m2_or_u2af1_vcf = select_first([U2AF1_wf.merged_vcf, Mutect2.filtered_vcf])
            File m2_or_u2af1_vcf_idx = select_first([U2AF1_wf.merged_vcf_idx, Mutect2.filtered_vcf_idx])
        }
    }

    Array[File] m2_vcfs_out = select_first([m2_or_u2af1_vcf, []])
    Array[File] m2_vcfs_idx_out = select_first([m2_or_u2af1_vcf_idx, []])
    Array[File] pon_input_vcfs = flatten([vcf_idx_pairs[0], m2_vcfs_out])
    Array[File] pon_input_vcfs_idx = flatten([vcf_idx_pairs[1], m2_vcfs_idx_out])

    if (!m2_only) {
        call m2.SplitIntervals {
            input:
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                ref_dict = ref_dict,
                scatter_count = select_first([create_panel_scatter_count, 24]),
                split_intervals_extra_args = "--dont-mix-contigs --min-contig-size " + contig_size,
                runtime_params = standard_runtime
        }

        scatter (subintervals in SplitIntervals.interval_files ) {
            call CreatePanel {
                input:
                    input_vcfs = pon_input_vcfs,
                    input_vcfs_idx = pon_input_vcfs_idx,
                    intervals = subintervals,
                    ref_fasta = ref_fasta,
                    ref_fai = ref_fai,
                    ref_dict = ref_dict,
                    gnomad = gnomad,
                    gnomad_idx = gnomad_idx,
                    output_vcf_name = pon_name,
                    create_pon_extra_args = create_pon_extra_args,
                    mem_mb = pon_mem,
                    mem_pad = command_mem_padding,
                    runtime_params = standard_runtime
            }
        }

        call m2.MergeVCFs {
            input:
                input_vcfs = CreatePanel.output_vcf,
                input_vcf_indices = CreatePanel.output_vcf_index,
                output_name = pon_name,
                compress = select_first([compress, false]),
                runtime_params = standard_runtime
        }
    }

    output {
        File? pon = MergeVCFs.merged_vcf
        File? pon_idx = MergeVCFs.merged_vcf_idx
        Array[File] normal_calls = pon_input_vcfs
        Array[File] normal_calls_idx = pon_input_vcfs_idx
    }
}

task CreatePanel {
    input {
        File intervals
        Array[File] input_vcfs
        Array[File] input_vcfs_idx
        File ref_fasta
        File ref_fai
        File ref_dict
        String output_vcf_name
        File gnomad
        File gnomad_idx
        String? create_pon_extra_args
        Int mem_mb = 5000
        Int mem_pad = 1000

        # runtime
        Runtime runtime_params
    }

    Int machine_mem = mem_mb
    Int cpu_mult = if runtime_params.cpu > 1 then runtime_params.cpu - 1 else 1
    Int command_mem = machine_mem - mem_pad

        parameter_meta{
            gnomad: {localization_optional: true}
            gnomad_idx: {localization_optional: true}
        }

    command {
        set -e
        export GATK_LOCAL_JAR=~{default="/gatk/gatk.jar" runtime_params.gatk_override}

        gatk --java-options "-Xmx~{command_mem}m -Xms~{command_mem - 1000}m" GenomicsDBImport --genomicsdb-workspace-path pon_db -R ~{ref_fasta} -V ~{sep=' -V ' input_vcfs} -L ~{intervals}

        gatk --java-options "-Xmx~{command_mem}m -Xms~{command_mem - 1000}m" CreateSomaticPanelOfNormals -R ~{ref_fasta} --germline-resource ~{gnomad} \
            -V gendb://pon_db -O ~{output_vcf_name}.vcf ~{create_pon_extra_args}
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
        File output_vcf = "~{output_vcf_name}.vcf"
        File output_vcf_index = "~{output_vcf_name}.vcf.idx"
    }
}
