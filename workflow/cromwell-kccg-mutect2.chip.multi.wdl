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
        # whitelist settings
        Int whitelist_mem_mb = 10000
        Int whitelist_disk_space = 300
        Int whitelist_cpu = 1
        Boolean treat_missing_as_rare = true
        Boolean whitelist_genome = true
        Boolean whitelist_use_ensembl_annotation = false
        String gnomad_pop = "AF"
        String whitelist_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/whitelist_filter@sha256:1f1f83f8241f40fbd1f21b19e2ccbdc184984fd9ec0b0a7bdfa97b8a73fed8a4"  # :latest
        File whitelist_archive
        File ref_fasta
        # common settings
        String ref_name = "hg38"
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
        Boolean use_sys_tmp_dir = true
    }

    Array[Array[String]] input_vcfs = read_tsv(input_vcf_list)

    scatter( row in input_vcfs ) {
        String tumor_sample_name = row[0]
        File input_vcf = row[1]

        call m2c.Mutect2CHIP_CHIP {
            input:
                input_vcf = input_vcf,
                tumor_sample_name = tumor_sample_name,
                annovar_mem_mb = annovar_mem_mb,
                annovar_disk_space = annovar_disk_space,
                annovar_cpu = annovar_cpu,
                annovar_docker = annovar_docker,
                annovar_archive = annovar_archive,
                whitelist_mem_mb = whitelist_mem_mb,
                whitelist_disk_space = whitelist_disk_space,
                whitelist_cpu = whitelist_cpu,
                treat_missing_as_rare = treat_missing_as_rare,
                whitelist_genome = whitelist_genome,
                whitelist_use_ensembl_annotation = whitelist_use_ensembl_annotation,
                gnomad_pop = gnomad_pop,
                whitelist_docker = whitelist_docker,
                whitelist_archive = whitelist_archive,
                ref_fasta = ref_fasta,
                ref_name = ref_name,
                preemptible = preemptible,
                max_retries = max_retries,
                small_task_cpu = small_task_cpu,
                small_task_mem = small_task_mem,
                small_task_disk = small_task_disk,
                command_mem_padding = command_mem_padding,
                boot_disk_size = boot_disk_size,
                emergency_extra_disk = emergency_extra_disk,
                use_sys_tmp_dir = use_sys_tmp_dir
        }
    }

    output {
        Array[File] out_whitelist_annovar_vcf = Mutect2CHIP_CHIP.out_whitelist_annovar_vcf
        Array[File] out_whitelist_annovar_table = Mutect2CHIP_CHIP.out_whitelist_annovar_table
        Array[File?] out_whitelist_annovar_output_refgene_variant_function = Mutect2CHIP_CHIP.out_whitelist_annovar_output_refgene_variant_function
        Array[File?] out_whitelist_annovar_output_ensgene_variant_function = Mutect2CHIP_CHIP.out_whitelist_annovar_output_ensgene_variant_function
        Array[File?] out_whitelist_annovar_output_refgene_variant_exonic_function = Mutect2CHIP_CHIP.out_whitelist_annovar_output_refgene_variant_exonic_function
        Array[File?] out_whitelist_annovar_output_ensgene_variant_exonic_function = Mutect2CHIP_CHIP.out_whitelist_annovar_output_ensgene_variant_exonic_function
        Array[File] out_whitelist_filter_output_allvariants_csv = Mutect2CHIP_CHIP.out_whitelist_filter_output_allvariants_csv
        Array[File] out_whitelist_filter_output_allvariantsfiltered_csv = Mutect2CHIP_CHIP.out_whitelist_filter_output_allvariantsfiltered_csv
        Array[File] out_whitelist_filter_output_exonicsplicingvariants_csv = Mutect2CHIP_CHIP.out_whitelist_filter_output_exonicsplicingvariants_csv
        Array[File] out_whitelist_filter_output_chiptranscriptvariants_csv = Mutect2CHIP_CHIP.out_whitelist_filter_output_chiptranscriptvariants_csv
        Array[File] out_whitelist_filter_output_chiptranscriptvariantsfiltered_csv = Mutect2CHIP_CHIP.out_whitelist_filter_output_chiptranscriptvariantsfiltered_csv
        Array[File] out_whitelist_filter_output_putativefilter_csv = Mutect2CHIP_CHIP.out_whitelist_filter_output_putativefilter_csv
    }
}