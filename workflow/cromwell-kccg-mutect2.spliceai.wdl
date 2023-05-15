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
# workflow developed by Alex Bick.                                #
#                                                                 #
# This pipeline has been developed for use by the Kinghorn        #
# Centre for Clinical Genomics and the Garvan Institute for       # 
# Medical Research.                                               #
#                                                                 #
# Author: Michael Geaghan (micgea)                                #
# Created: 2021/08/13                                             #
# =============================================================== #

import "cromwell-kccg-mutect2.wdl" as m2

workflow Mutect2CHIP_SpliceAI {
    input {
        File input_vcf
        File ref_fasta
        File? spliceai_annotation_file
        String spliceai_annotation_string = 'grch38'
        Int spliceai_max_dist = 50
        Boolean spliceai_mask = false
        String spliceai_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/spliceai@sha256:b3dfb27959e8a8ef6a7d0bca1562f86765ba7d4dffd691b83aa94cf733785a8d"  # :v1.3
        Int spliceai_disk = 100
        Int spliceai_mem_mb = 16000
        Int spliceai_cpu = 4
        # runtime parameters
        Int? preemptible
        Int? max_retries
        Int small_task_cpu = 4
        Int small_task_mem = 4000
        Int small_task_disk = 100
        Int command_mem_padding = 1000
        Int boot_disk_size = 12
        # Use as a last resort to increase the disk given to every task in case of ill behaving data
        Int? emergency_extra_disk
    }

    Int preemptible_or_default = select_first([preemptible, 2])
    Int max_retries_or_default = select_first([max_retries, 2])

    # This is added to every task as padding, should increase if systematically you need more disk for every call
    Int disk_pad = 10 + select_first([emergency_extra_disk,0])

    Runtime standard_runtime = {
        "gatk_docker": "",
        "max_retries": max_retries_or_default,
        "preemptible": preemptible_or_default,
        "cpu": small_task_cpu,
        "machine_mem": small_task_mem,
        "command_mem": small_task_mem - command_mem_padding,
        "disk": small_task_disk + disk_pad,
        "boot_disk_size": boot_disk_size
    }

    call m2.SpliceAI {
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
        File spliceai_output_vcf = SpliceAI.spliceai_output_vcf
    }
}