version 1.0

# =============================================================== #
# cromwell-kccg-mutect2-chip                                      #
#                                                                 #
# This workflow detects mutations associated with CHIP            #
# (clonal haematopoiesis of indeterminate potential) using the    #
# GATK4 Mutect2 somatic variant calling pipeline and the          #
# Cromwell workflow engine.                                       #
#                                                                 #
# This script will call the annovar annotation pipeline on        #
# multiple samples given an input TSV file specifying the input   #
# vcf files.                                                      #
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

import "cromwell-kccg-mutect2.annovar.wdl" as m2a

workflow Mutect2CHIP_Annovar_Multi {
    input {
        # input vcf
        File input_vcf_list
        # annovar settings
        Int annovar_mem_mb = 4000
        Int annovar_disk = 100
        Int annovar_tmp_disk = 200
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
    }

    Array[Array[String]] input_vcfs = read_tsv(input_vcf_list)

    scatter( row in input_vcfs ) {
        File input_vcf = row[0]

        call m2a.Mutect2CHIP_Annovar {
            input:
                input_vcf = input_vcf,
                annovar_mem_mb = annovar_mem_mb,
                annovar_disk = annovar_disk,
                annovar_tmp_disk = annovar_tmp_disk,
                annovar_cpu = annovar_cpu,
                annovar_docker = annovar_docker,
                annovar_archive = annovar_archive,
                ref_name = ref_name,
                annovar_protocols = annovar_protocols,
                annovar_operations = annovar_operations,
                preemptible = preemptible,
                max_retries = max_retries,
                small_task_cpu = small_task_cpu,
                small_task_mem = small_task_mem,
                small_task_disk = small_task_disk,
                command_mem_padding = command_mem_padding,
                boot_disk_size = boot_disk_size,
                emergency_extra_disk = emergency_extra_disk,
        }
    }

    output {
        Array[File] out_annovar_vcf = Mutect2CHIP_Annovar.out_annovar_vcf 
        Array[File] out_annovar_table = Mutect2CHIP_Annovar.out_annovar_table 
    }
}