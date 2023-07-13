version 1.0

# =============================================================== #
# chip.wdl                                                        #
#                                                                 #
# This workflow detects mutations associated with CHIP            #
# (clonal haematopoiesis of indeterminate potential).             #
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
# Created: 2023/04/04                                             #
# =============================================================== #

import "annovar.wdl" as Annovar

workflow CHIP {
    input {
        # input vcf
        File input_vcf
        String tumor_sample_name
        # annovar settings
        Int annovar_mem_mb = 4000
        Int annovar_disk = 100
        Int annovar_tmp_disk = 200
        Int annovar_cpu = 1
        String annovar_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/perl@sha256:1f35086e2ff48dace3b3edeaa2ad1faf1e44c0612e00f00ea0fc1830b576a261"  # :5.34.0
        File annovar_archive
        # whitelist settings
        Int whitelist_mem_mb = 10000
        Int whitelist_disk = 300
        Int whitelist_cpu = 1
        Boolean treat_missing_as_rare = true
        Boolean whitelist_genome = true
        Boolean whitelist_use_ensembl_annotation = false
        String gnomad_pop = "AF"
        String whitelist_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/whitelist_filter@sha256:9cd77186c23a0b256a0928c5a4087b8378c234cb0754f35557cf9ec6d4aa544d"  # :latest
        File whitelist_archive
        File ref_fasta
        # common settings
        String ref_name = "hg38"
        # runtime parameters
        Int? preemptible
        Int? max_retries
        Int boot_disk_size = 12
    }

    Int preemptible_or_default = select_first([preemptible, 2])
    Int max_retries_or_default = select_first([max_retries, 2])

    AnnotationRuntime standard_runtime = {
        "max_retries": max_retries_or_default,
        "preemptible": preemptible_or_default,
        "boot_disk_size": boot_disk_size
    }

    String sample_id = basename(basename(input_vcf, ".gz"), ".vcf")
    String whitelist_annovar_protocols = if (whitelist_genome) then "refGene,ensGene,gnomad211_genome,gnomad211_exome" else "refGene,ensGene,gnomad211_exome"
    String whitelist_annovar_operations = if (whitelist_genome) then "g,g,f,f" else "g,g,f"
    String additional_arguments = '-argument "-exonicsplicing -transcript_function -separate,-exonicsplicing -transcript_function -separate,'
    String additional_arguments_final = if (whitelist_genome) then additional_arguments + ',"' else additional_arguments + '"'
    call Annovar.Annovar_task as WhitelistAnnovar {
        input:
            vcf_input = input_vcf,
            sample_id = sample_id,
            mem_mb = annovar_mem_mb,
            annovar_disk_space = annovar_disk,
            annovar_tmp_disk_space = annovar_tmp_disk,
            cpu = annovar_cpu,
            annovar_docker = annovar_docker,
            annovar_archive = annovar_archive,
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

    call WhitelistFilter {
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
            whitelist_filter_archive = whitelist_archive,
            runtime_params = standard_runtime
    }

    output {
        File out_whitelist_annovar_vcf = WhitelistAnnovar.annovar_output_file_vcf
        File out_whitelist_annovar_table = WhitelistAnnovar.annovar_output_file_table
        File? out_whitelist_annovar_output_refgene_variant_function = WhitelistAnnovar.annovar_output_refgene_variant_function
        File? out_whitelist_annovar_output_ensgene_variant_function = WhitelistAnnovar.annovar_output_ensgene_variant_function
        File? out_whitelist_annovar_output_refgene_variant_exonic_function = WhitelistAnnovar.annovar_output_refgene_variant_exonic_function
        File? out_whitelist_annovar_output_ensgene_variant_exonic_function = WhitelistAnnovar.annovar_output_ensgene_variant_exonic_function
        File out_whitelist_filter_output_csv = WhitelistFilter.whitelist_filter_output_csv
        File out_whitelist_filter_output_allvariants_csv = WhitelistFilter.whitelist_filter_output_allvariants_csv
        File out_whitelist_filter_output_allvariantsfiltered_csv = WhitelistFilter.whitelist_filter_output_allvariantsfiltered_csv
        File out_whitelist_filter_output_exonicsplicingvariants_csv = WhitelistFilter.whitelist_filter_output_exonicsplicingvariants_csv
        File out_whitelist_filter_output_chiptranscriptvariants_csv = WhitelistFilter.whitelist_filter_output_chiptranscriptvariants_csv
        File out_whitelist_filter_output_chiptranscriptvariantsfiltered_csv = WhitelistFilter.whitelist_filter_output_chiptranscriptvariantsfiltered_csv
        File out_whitelist_filter_output_putativefilter_csv = WhitelistFilter.whitelist_filter_output_putativefilter_csv
    }
}

task WhitelistFilter {
    input {
      Int mem_mb = 10000
      Int whitelist_filter_disk_space = 300
      Int cpu = 1
      String whitelist_filter_docker
      String tumor_sample_name
      Boolean use_ensembl_annotation = false
      String gnomad_source = "genome,exome"
      String gnomad_pop = "AF"
      Boolean treat_missing_as_rare = true
      File txt_input
      File vcf_input
      File var_func_input
      File var_exonic_func_input
      String ref_name
      File ref_fasta
      File whitelist_filter_archive
      AnnotationRuntime runtime_params
    }

    String file_prefix = basename(txt_input, "_multianno.txt")
    String treat_missing_as_rare_str = if (treat_missing_as_rare) then "TRUE" else "FALSE"
    String ensembl_refseq = if (use_ensembl_annotation) then "ensembl" else "refseq"

    command {
      set -euo pipefail

      tar -xzvf ~{whitelist_filter_archive}

      SCRIPT_DIR=$PWD
      cd whitelist

      Rscript ./whitelist_filter_rscript.R \
        ~{vcf_input} \
        ~{txt_input} \
        ~{var_func_input} \
        ~{var_exonic_func_input} \
        ~{ensembl_refseq} \
        ~{tumor_sample_name} \
        ~{gnomad_source} \
        ~{gnomad_pop} \
        ~{treat_missing_as_rare_str} \
        ./chip_variant_definitions.csv \
        ~{ref_fasta} \
        ./somaticism_filter_transcripts.txt

      mv *_variants.csv *_variants.*.csv $SCRIPT_DIR/
      cd $SCRIPT_DIR
      mkdir -p chip_csv_outputs
      mv *_variants.csv *_variants.*.csv chip_csv_outputs/
      cp chip_csv_outputs/*.putative_filter.csv ./~{file_prefix}.chip.csv
    }

    runtime {
      docker: whitelist_filter_docker
      bootDiskSizeGb: runtime_params.boot_disk_size
      memory: mem_mb + " MB"
      disks: "local-disk " + whitelist_filter_disk_space + " HDD"
      preemptible: runtime_params.preemptible
      maxRetries: runtime_params.max_retries
      cpu: cpu
    }

    output {
      File whitelist_filter_output_csv = file_prefix + ".chip.csv"
      File whitelist_filter_output_allvariants_csv = "chip_csv_outputs/" + file_prefix + ".all_variants.csv"
      File whitelist_filter_output_allvariantsfiltered_csv = "chip_csv_outputs/" + file_prefix + ".all_variants.pre_filtered.csv"
      File whitelist_filter_output_exonicsplicingvariants_csv = "chip_csv_outputs/" + file_prefix + ".exonic_splicing_variants.csv"
      File whitelist_filter_output_chiptranscriptvariants_csv = "chip_csv_outputs/" + file_prefix + ".chip_transcript_variants.csv"
      File whitelist_filter_output_chiptranscriptvariantsfiltered_csv = "chip_csv_outputs/" + file_prefix + ".chip_transcript_variants.filtered.csv"
      File whitelist_filter_output_putativefilter_csv = "chip_csv_outputs/" + file_prefix + ".chip_transcript_variants.filtered.putative_filter.csv"
    }
}
