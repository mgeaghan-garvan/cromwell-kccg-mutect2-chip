version 1.0

# =============================================================== #
# cromwell-kccg-mutect2-chip                                      #
#                                                                 #
# This workflow detects mutations associated with CHIP            #
# (clonal haematopoiesis of indeterminate potential) using the    #
# GATK4 Mutect2 somatic variant calling pipeline and the          #
# Cromwell workflow engine.                                       #
#                                                                 #
# This script will call the CHIP cromwell pipeline on multiple    #
# samples given an input TSV file specifying the input VCF        #
# files.                                                          #
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

import "cromwell-kccg-mutect2.chip.wdl" as m2c

workflow Mutect2CHIP_CHIP_Multi {
    input {
        # input vcf
        File input_vcf_list
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
        Int small_task_cpu = 2
        Int small_task_mem = 4000
        Int small_task_disk = 100
        Int boot_disk_size = 12
        # Use as a last resort to increase the disk given to every task in case of ill behaving data
        Int? emergency_extra_disk
    }

    # TODO: test this - I'm not sure if this works for single-column files. It might only yield Array[String]
    Array[Array[String]] input_vcfs = read_tsv(input_vcf_list)

    scatter( row in input_vcfs ) {
        File input_vcf = row[0]

        call m2c.Mutect2CHIP_CHIP {
            input:
                input_vcf = input_vcf,
                annovar_mem_mb = annovar_mem_mb,
                annovar_disk_space = annovar_disk_space,
                annovar_cpu = annovar_cpu,
                annovar_docker = annovar_docker,
                annovar_archive = annovar_archive,
                annovar_protocols = annovar_protocols,
                annovar_operations = annovar_operations,
                whitelist_mem_mb = whitelist_mem_mb,
                whitelist_disk_space = whitelist_disk_space,
                whitelist_cpu = whitelist_cpu,
                whitelist_docker = whitelist_docker,
                whitelist_archive = whitelist_archive,
                ref_name = ref_name,
                gatk_docker = gatk_docker,
                gatk_override = gatk_override,
                preemptible = preemptible,
                max_retries = max_retries,
                small_task_cpu = small_task_cpu,
                small_task_mem = small_task_mem,
                small_task_disk = small_task_disk,
                boot_disk_size = boot_disk_size,
                emergency_extra_disk = emergency_extra_disk
        }
    }

    output {
        Array[File] out_annovar_vcf = Mutect2CHIP_CHIP.out_annovar_vcf 
        Array[File] out_annovar_table = Mutect2CHIP_CHIP.out_annovar_table 
        Array[File?] out_whitelist_count = Mutect2CHIP_CHIP.out_whitelist_count 
        Array[File?] out_whitelist_all_variants = Mutect2CHIP_CHIP.out_whitelist_all_variants 
        Array[File?] out_whitelist = Mutect2CHIP_CHIP.out_whitelist 
        Array[File?] out_whitelist_manual_review = Mutect2CHIP_CHIP.out_whitelist_manual_review 
    }
}