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

workflow Mutect2CHIP_Annovar {
    input {
        # input vcf
        File input_vcf
        # annovar settings
        Int annovar_mem_mb = 4000
        Int annovar_disk_space = 300
        Int annovar_cpu = 1
        String annovar_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/perl@sha256:1f35086e2ff48dace3b3edeaa2ad1faf1e44c0612e00f00ea0fc1830b576a261"  # :5.34.0
        File annovar_archive
        String ref_name = "hg38"
        String annovar_protocols = "cosmic70"
        String annovar_operations = "f"
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
        Boolean use_tmp_dir = true
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

    String sample_id = basename(basename(input_vcf, ".gz"), ".vcf")
    call m2.Annovar {
        input:
            mem_mb = annovar_mem_mb,
            annovar_disk_space = annovar_disk_space,
            cpu = annovar_cpu,
            annovar_docker = annovar_docker,
            sample_id = sample_id,
            vcf_input = input_vcf,
            annovar_archive = annovar_archive,
            ref_name = ref_name,
            annovar_protocols = annovar_protocols,
            annovar_operations = annovar_operations,
            use_tmp_dir = use_tmp_dir,
            runtime_params = standard_runtime
    }

    output {
        File out_annovar_vcf = Annovar.annovar_output_file_vcf
        File out_annovar_table = Annovar.annovar_output_file_table
    }
}