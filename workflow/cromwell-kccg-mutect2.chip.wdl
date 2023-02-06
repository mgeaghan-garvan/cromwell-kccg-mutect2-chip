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
        String tumor_sample_name
        # annovar settings
        Int annovar_cpu = 1
        String annovar_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/perl@sha256:1f35086e2ff48dace3b3edeaa2ad1faf1e44c0612e00f00ea0fc1830b576a261"  # :5.34.0
        File annovar_archive
        # whitelist settings
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
        Int annovar_mem_mb = 4000
        Int annovar_disk = 100
        Int annovar_tmp_disk = 200
        Int whitelist_mem_mb = 10000
        Int whitelist_disk = 300
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

    File annovar_archive_file = select_first([annovar_archive, "ANNOVAR_ARCHIVE_NOT_SUPPLIED"])
    File whitelist_archive_file = select_first([whitelist_archive, "WHITELIST_FILTER_ARCHIVE_NOT_SUPPLIED"])
    String sample_id = basename(basename(input_vcf, ".gz"), ".vcf")
    String whitelist_annovar_protocols = if (whitelist_genome) then "refGene,ensGene,gnomad211_genome,gnomad211_exome" else "refGene,ensGene,gnomad211_exome"
    String whitelist_annovar_operations = if (whitelist_genome) then "g,g,f,f" else "g,g,f"
    String additional_arguments = '-argument "-exonicsplicing -transcript_function -separate,-exonicsplicing -transcript_function -separate,'
    String additional_arguments_final = if (whitelist_genome) then additional_arguments + ',"' else additional_arguments + '"'
    call m2.Annovar as WhitelistAnnovar {
        input:
            mem_mb = annovar_mem_mb,
            annovar_disk_space = annovar_disk,
            annovar_tmp_disk_space = annovar_tmp_disk,
            cpu = annovar_cpu,
            annovar_docker = annovar_docker,
            sample_id = sample_id,
            vcf_input = input_vcf,
            annovar_archive = annovar_archive_file,
            ref_name = ref_name,
            label = "whitelist_annovar_out",
            annovar_protocols = whitelist_annovar_protocols,
            annovar_operations = whitelist_annovar_operations,
            annovar_additional_arguments = additional_arguments_final,
            runtime_params = standard_runtime
    }

    String whitelist_exome_only_or_both = if (whitelist_genome) then "genome,exome" else "exome"
    File var_func_input = if (whitelist_use_ensembl_annotation) then select_first([WhitelistAnnovar.annovar_output_ensgene_variant_function]) else select_first([WhitelistAnnovar.annovar_output_refgene_variant_function])
    File var_exonic_func_input = if (whitelist_use_ensembl_annotation) then select_first([WhitelistAnnovar.annovar_output_ensgene_variant_exonic_function]) else select_first([WhitelistAnnovar.annovar_output_refgene_variant_exonic_function])

    call m2.WhitelistFilter {
        input:
            mem_mb = whitelist_mem_mb,
            whitelist_filter_disk_space = whitelist_disk,
            cpu = whitelist_cpu,
            whitelist_filter_docker = whitelist_docker,
            tumor_sample_name = tumor_sample_name,
            use_ensembl_annotation = whitelist_use_ensembl_annotation,
            gnomad_source = whitelist_exome_only_or_both,
            gnomad_pop = gnomad_pop,
            treat_missing_as_rare = treat_missing_as_rare,
            txt_input = WhitelistAnnovar.annovar_output_file_table,
            vcf_input = WhitelistAnnovar.annovar_output_file_vcf,
            var_func_input = var_func_input,
            var_exonic_func_input = var_exonic_func_input,
            ref_name = ref_name,
            ref_fasta = ref_fasta,
            whitelist_filter_archive = whitelist_archive_file,
            runtime_params = standard_runtime
    }

    output {
        File out_whitelist_annovar_vcf = WhitelistAnnovar.annovar_output_file_vcf
        File out_whitelist_annovar_table = WhitelistAnnovar.annovar_output_file_table
        File? out_whitelist_annovar_output_refgene_variant_function = WhitelistAnnovar.annovar_output_refgene_variant_function
        File? out_whitelist_annovar_output_ensgene_variant_function = WhitelistAnnovar.annovar_output_ensgene_variant_function
        File? out_whitelist_annovar_output_refgene_variant_exonic_function = WhitelistAnnovar.annovar_output_refgene_variant_exonic_function
        File? out_whitelist_annovar_output_ensgene_variant_exonic_function = WhitelistAnnovar.annovar_output_ensgene_variant_exonic_function
        File out_whitelist_filter_output_allvariants_csv = WhitelistFilter.whitelist_filter_output_allvariants_csv
        File out_whitelist_filter_output_allvariantsfiltered_csv = WhitelistFilter.whitelist_filter_output_allvariantsfiltered_csv
        File out_whitelist_filter_output_exonicsplicingvariants_csv = WhitelistFilter.whitelist_filter_output_exonicsplicingvariants_csv
        File out_whitelist_filter_output_chiptranscriptvariants_csv = WhitelistFilter.whitelist_filter_output_chiptranscriptvariants_csv
        File out_whitelist_filter_output_chiptranscriptvariantsfiltered_csv = WhitelistFilter.whitelist_filter_output_chiptranscriptvariantsfiltered_csv
        File out_whitelist_filter_output_putativefilter_csv = WhitelistFilter.whitelist_filter_output_putativefilter_csv
    }
}