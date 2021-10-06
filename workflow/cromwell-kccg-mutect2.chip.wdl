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
# The CHIP detection stage is adapted from the                    #
# Annovar Whitelist Filter developed by Charlie Condon.           #
#                                                                 #
# This pipeline has been developed for use by the Kinghorn        #
# Centre for Clinical Genomics and the Garvan Institute for       # 
# Medical Research.                                               #
#                                                                 #
# Author: Michael Geaghan (micgea)                                #
# Created: 2021/08/13                                             #
# =============================================================== #

import "cromwell-kccg-mutect2.wdl" as m2

workflow Mutect2CHIP_CHIP {
    input {
        # input vcf
        File input_vcf
        # annovar settings
        Int annovar_mem_mb = 4000
        Int annovar_disk_space = 300
        Int annovar_cpu = 1
        String annovar_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/perl@sha256:1f35086e2ff48dace3b3edeaa2ad1faf1e44c0612e00f00ea0fc1830b576a261"  # :5.34.0
        File annovar_archive
        String annovar_protocols = "refGene,cosmic70"
        String annovar_operations = "g,f"
        # whitelist settings
        Int whitelist_mem_mb = 10000
        Int whitelist_disk_space = 300
        Int whitelist_cpu = 1
        String whitelist_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/whitelist_filter@sha256:3e3868fbb7e58e6f9550cf15c046e6c004a28b8e98b1008224a272d82a4dc357"  # :latest
        File whitelist_archive
        # common settings
        String ref_name = "hg38"
        # runtime parameters
        String gatk_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/gatk@sha256:0359ae4f32f2f541ca86a8cd30ef730bbaf8c306b9d53d2d520262d3e84b3b2b"  # :4.2.1.0
        File? gatk_override
        Int? preemptible
        Int? max_retries
        Int small_task_cpu = 4
        Int small_task_mem = 4000
        Int small_task_disk = 100
        Int boot_disk_size = 12
        # Use as a last resort to increase the disk given to every task in case of ill behaving data
        Int? emergency_extra_disk
    }

    Int preemptible_or_default = select_first([preemptible, 2])
    Int max_retries_or_default = select_first([max_retries, 2])

    # If no tar is provided, the task downloads one from broads ftp server
    Int gatk_override_size = if defined(gatk_override) then ceil(size(gatk_override, "GB")) else 0

    # This is added to every task as padding, should increase if systematically you need more disk for every call
    Int disk_pad = 10 + gatk_override_size + select_first([emergency_extra_disk,0])

    Runtime standard_runtime = {
        "gatk_docker": gatk_docker,
        "gatk_override": gatk_override,
        "max_retries": max_retries_or_default,
        "preemptible": preemptible_or_default,
        "cpu": small_task_cpu,
        "machine_mem": small_task_mem,
        "command_mem": small_task_mem - 500,
        "disk": small_task_disk + disk_pad,
        "boot_disk_size": boot_disk_size
    }

    File annovar_archive_file = select_first([annovar_archive, "ANNOVAR_ARCHIVE_NOT_SUPPLIED"])
    File whitelist_archive_file = select_first([whitelist_archive, "WHITELIST_FILTER_ARCHIVE_NOT_SUPPLIED"])
    String sample_id = basename(basename(input_vcf, ".gz"), ".vcf")
    call m2.Annovar {
        input:
            mem_mb = annovar_mem_mb,
            annovar_disk_space = annovar_disk_space,
            cpu = annovar_cpu,
            annovar_docker = annovar_docker,
            sample_id = sample_id,
            vcf_input = input_vcf,
            annovar_archive = annovar_archive_file,
            ref_name = ref_name,
            runtime_params = standard_runtime
    }

    call m2.WhitelistFilter {
        input:
            mem_mb = whitelist_mem_mb,
            whitelist_filter_disk_space = whitelist_disk_space,
            cpu = whitelist_cpu,
            whitelist_filter_docker = whitelist_docker,
            txt_input = Annovar.annovar_output_file_table,
            ref_name = ref_name,
            whitelist_filter_archive = whitelist_archive_file,
            runtime_params = standard_runtime
    }

    output {
        File out_annovar_vcf = Annovar.annovar_output_file_vcf
        File out_annovar_table = Annovar.annovar_output_file_table
        File? out_whitelist_count = WhitelistFilter.whitelist_filter_output_varcount_csv
        File? out_whitelist_all_variants = WhitelistFilter.whitelist_filter_output_allvariants_csv
        File? out_whitelist = WhitelistFilter.whitelist_filter_output_wl_csv
        File? out_whitelist_manual_review = WhitelistFilter.whitelist_filter_output_manual_review_csv
    }
}