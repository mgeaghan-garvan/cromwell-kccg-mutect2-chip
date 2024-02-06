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
        Int annovar_cpu = 1
        String annovar_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/annovar@sha256:842e9f88dd39999ee2129aeb992e8eced10ac2a33642d4b34d0f0c0254aa5035"  # :5.34.0
        Array[File] annovar_db_files
        String ref_name = "hg38"
        # chip annotation settings
        Int chip_mem_mb = 10000
        Int chip_disk = 300
        Int chip_cpu = 1
        Boolean treat_missing_as_rare = true
        Boolean use_gnomad_genome = true
        Boolean use_ensembl_annotation = false
        String gnomad_pop = "AF"
        String chip_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/whitelist_filter@sha256:9cd77186c23a0b256a0928c5a4087b8378c234cb0754f35557cf9ec6d4aa544d"  # :latest
        File ref_fasta
        File chip_mutations_csv
        File somaticism_filter_transcripts
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

    call CHIPPreFilter {
      input:
        input_vcf = input_vcf,
        chip_mutations_csv = chip_mutations_csv,
        ref_fasta = ref_fasta,
        mem_mb = chip_mem_mb,
        disk_space = chip_disk,
        cpu = chip_cpu,
        docker = chip_docker,
        runtime_params = standard_runtime
    }

    String chip_annovar_protocols = if (use_gnomad_genome) then "refGene,ensGene,gnomad211_genome,gnomad211_exome" else "refGene,ensGene,gnomad211_exome"
    String chip_annovar_operations = if (use_gnomad_genome) then "g,g,f,f" else "g,g,f"
    String additional_arguments = '-argument "-exonicsplicing -transcript_function -separate,-exonicsplicing -transcript_function -separate,'
    String additional_arguments_final = if (use_gnomad_genome) then additional_arguments + ',"' else additional_arguments + '"'
    call Annovar.Annovar_task as CHIPAnnovar {
        input:
            vcf_input = CHIPPreFilter.out_chip_genes_norm_no_info_filtered_vcf,
            sample_id = sample_id,
            mem_mb = annovar_mem_mb,
            annovar_disk_space = annovar_disk,
            cpu = annovar_cpu,
            annovar_docker = annovar_docker,
            annovar_db_files = annovar_db_files,
            ref_name = ref_name,
            label = "chip_annovar_out",
            annovar_protocols = chip_annovar_protocols,
            annovar_operations = chip_annovar_operations,
            annovar_additional_arguments = additional_arguments_final,
            runtime_params = standard_runtime
    }

    File var_func_input = if (use_ensembl_annotation) then select_first([CHIPAnnovar.annovar_output_ensgene_variant_function]) else select_first([CHIPAnnovar.annovar_output_refgene_variant_function])
    File var_exonic_func_input = if (use_ensembl_annotation) then select_first([CHIPAnnovar.annovar_output_ensgene_variant_exonic_function]) else select_first([CHIPAnnovar.annovar_output_refgene_variant_exonic_function])

    call CHIPAnnotation {
        input:
            annovar_txt_input = CHIPAnnovar.annovar_output_file_table,
            annovar_variant_function = var_func_input,
            annovar_exonic_function = var_exonic_func_input,
            chip_mutations_csv = chip_mutations_csv,
            somaticism_filter_transcripts = somaticism_filter_transcripts,
            seq_context_tsv = CHIPPreFilter.out_seq_tsv,
            input_vcf_header = CHIPPreFilter.out_chip_genes_norm_no_info_filtered_vcf_header,
            sample_id = tumor_sample_name,
            treat_missing_as_rare = treat_missing_as_rare,
            use_gnomad_genome = use_gnomad_genome,
            use_ensembl_annotation = use_ensembl_annotation,
            gnomad_pop = gnomad_pop,
            mem_mb = chip_mem_mb,
            disk_space = chip_disk,
            cpu = chip_cpu,
            docker = chip_docker,
            runtime_params = standard_runtime
    }

    call FinaliseCHIPFilter {
        input:
            input_vcf = CHIPPreFilter.gzipped_vcf,
            input_vcf_idx = CHIPPreFilter.gzipped_vcf_idx,
            non_chip_gene_annotations_vcf = CHIPPreFilter.out_non_chip_genes_so_filtered_vcf_gz,
            non_chip_gene_annotations_vf_idx = CHIPPreFilter.out_non_chip_genes_so_filtered_vcf_gz_idx,
            chip_gene_annotations_vcf = CHIPAnnotation.chip_vcf,
            chip_gene_annotations_vcf_idx = CHIPAnnotation.chip_vcf_idx,
            sample_id = tumor_sample_name,
            mem_mb = chip_mem_mb,
            disk_space = chip_disk,
            cpu = chip_cpu,
            docker = chip_docker,
            runtime_params = standard_runtime
    }

    output {
        File out_vcf = FinaliseCHIPFilter.chip_vcf
        File out_vcf_idx = FinaliseCHIPFilter.chip_vcf_idx
        File chip_vcf = CHIPAnnotation.chip_vcf
        File chip_vcf_idx = CHIPAnnotation.chip_vcf_idx
        File chip_csv = CHIPAnnotation.chip_csv
        File chip_rdata = CHIPAnnotation.chip_rdata
    }
}

task CHIPPreFilter {
    input {
      File input_vcf
      File chip_mutations_csv
      File ref_fasta
      Int mem_mb = 8000
      Int disk_space = 50
      Int cpu = 1
      String docker
      AnnotationRuntime runtime_params
    }

    String sample_id = basename(basename(input_vcf, ".gz"), ".vcf")
    String input_vcf_gz = sample_id + ".input_vcf.gz"
    String input_vcf_gz_idx = input_vcf_gz + ".tbi"
    String chip_genes_bed = "chip_genes.bed"
    String chip_genes_merged_bed = "chip_genes.merged.bed"
    String chip_genes_vcf = sample_id + ".chip_genes.vcf"
    String chip_genes_norm_vcf = sample_id + ".chip_genes.norm.vcf"
    String chip_genes_norm_no_info_vcf = sample_id + ".chip_genes.norm.no_info.vcf"
    String chip_genes_norm_no_info_filtered_vcf = sample_id + ".chip_genes.norm.no_info.filtered.vcf"
    String chip_genes_norm_no_info_filtered_vcf_header = chip_genes_norm_no_info_filtered_vcf + ".header"
    String chip_genes_norm_bed = sample_id + ".chip_genes.norm.bed"
    String non_chip_genes_so_filtered_vcf = sample_id + ".non_chip_genes.so.vcf"
    String non_chip_genes_so_filtered_vcf_gz = non_chip_genes_so_filtered_vcf + ".gz"
    String seq_tsv = sample_id + ".seq.tsv"

    command <<<
      set -euo pipefail

      # --- Step 0: Compress VCF if necessary and index ---
      VCF_SUFFIX=$(basename ~{input_vcf} | rev | cut -d "." -f 1 | rev)
      if [ "${VCF_SUFFIX}" != "gz" ]
      then
        bgzip -c ~{input_vcf} > ~{input_vcf_gz}
        tabix -s 1 -b 2 -e 2 ~{input_vcf_gz}
      else
        cp ~{input_vcf} ~{input_vcf_gz}
        bcftools index ~{input_vcf_gz}
      fi

      # --- Step 1: Turn the CHIP mutations CSV into a BED file of CHIP gene regions ---
      # This uses a quick R heredoc
      R --vanilla <<EOF
      library(tidyverse)
      df <- read_csv("~{chip_mutations_csv}")
      df <- df %>%
        select(c(chr, gene_genomic_start, gene_genomic_end)) %>%
        mutate(gene_genomic_start = gene_genomic_start - 1) %>%
        distinct %>%
        arrange(chr, gene_genomic_start, gene_genomic_end)
        write_tsv(df, "~{chip_genes_bed}", col_names = FALSE)
      EOF

      # Merge overlapping regions
      bedtools merge -i ~{chip_genes_bed} > ~{chip_genes_merged_bed}

      # --- Step 2: Filter the VCF to separate CHIP and non-CHIP genes ---
      # CHIP genes
      bedtools intersect -a ~{input_vcf_gz} -b ~{chip_genes_merged_bed} -wa -header > ~{chip_genes_vcf}
      # Normalise to split multi-allelic variants
      bcftools norm -m -any -o ~{chip_genes_norm_vcf} ~{chip_genes_vcf}
      # Strip INFO field
      bcftools annotate -x INFO -o ~{chip_genes_norm_no_info_vcf} ~{chip_genes_norm_vcf}
      
      # Non-CHIP genes -> sites-only VCF
      bedtools intersect -a ~{input_vcf_gz} -b ~{chip_genes_merged_bed} -wa -header -v | \
        cut -f 1-8 | \
        awk -v FS="\t" -v OFS="\t" '
          $0 ~ /^#/ { print $0 }
          $0 !~ /^#/ { print $1, $2, $3, $4, $5, $6, "non_chip_gene", "." }
        ' > ~{non_chip_genes_so_filtered_vcf}
      # Add new FILTER header
      NEW_FILTER_HEADER="##FILTER=<ID=non_chip_gene,Description=\"Variant is not in a CHIP gene region\">"
      LAST_FILTER_LINE=$(grep -n "^##FILTER" ~{non_chip_genes_so_filtered_vcf} | tail -n 1 | cut -d ":" -f 1)
      sed -i "${LAST_FILTER_LINE}a ${NEW_FILTER_HEADER}" ~{non_chip_genes_so_filtered_vcf}
      bgzip -c ~{non_chip_genes_so_filtered_vcf} > ~{non_chip_genes_so_filtered_vcf_gz}
      tabix -s 1 -b 2 -e 2 ~{non_chip_genes_so_filtered_vcf_gz}

      # --- Step 3: Create sequence context TSV file ---
      awk -v OFS="\t" '$0 !~ /^#/ { print $1, $2 - 11, $2 + length($4) + 9, $1":"$2":"$4":"$5 }' ~{chip_genes_norm_no_info_vcf} > ~{chip_genes_norm_bed}
      bedtools getfasta -fi ~{ref_fasta} -bed ~{chip_genes_norm_bed} -tab -nameOnly > ~{seq_tsv}

      # --- Step 4: Apply filters to CHIP genes VCF ---
      bcftools filter \
        -s chip_ad_filter_fail \
        -m + \
        -i 'FORMAT/AD[0:1]>=3' \
        ~{chip_genes_norm_no_info_vcf} | \
      bcftools filter \
        -s chip_dp_filter_fail \
        -m + \
        -i 'FORMAT/DP[0]>=20' | \
      bcftools filter \
        -s chip_f1r2_filter_fail \
        -m + \
        -e 'FORMAT/F1R2[0:*]<1' | \
      bcftools filter \
        -s chip_f2r1_filter_fail \
        -m + \
        -e 'FORMAT/F2R1[0:*]<1' > ~{chip_genes_norm_no_info_filtered_vcf}

      # --- Step 5: Get VCF header ---
      grep "^#" ~{chip_genes_norm_no_info_filtered_vcf} > ~{chip_genes_norm_no_info_filtered_vcf_header}
    >>>

    runtime {
      docker: docker
      bootDiskSizeGb: runtime_params.boot_disk_size
      memory: mem_mb + " MB"
      disks: "local-disk " + disk_space + " HDD"
      preemptible: runtime_params.preemptible
      maxRetries: runtime_params.max_retries
      cpu: cpu
    }

    output {
      File gzipped_vcf = input_vcf_gz
      File gzipped_vcf_idx = input_vcf_gz_idx
      File out_chip_genes_norm_no_info_filtered_vcf = chip_genes_norm_no_info_filtered_vcf
      File out_chip_genes_norm_no_info_filtered_vcf_header = chip_genes_norm_no_info_filtered_vcf_header
      File out_non_chip_genes_so_filtered_vcf_gz = non_chip_genes_so_filtered_vcf_gz
      File out_non_chip_genes_so_filtered_vcf_gz_idx = non_chip_genes_so_filtered_vcf_gz + ".tbi"
      File out_seq_tsv = seq_tsv
    }
}

task CHIPAnnotation {
  input {
      File annovar_txt_input
      File annovar_variant_function
      File annovar_exonic_function
      File chip_mutations_csv
      File somaticism_filter_transcripts
      File seq_context_tsv
      File input_vcf_header
      String sample_id
      Boolean treat_missing_as_rare = true
      Boolean use_gnomad_genome = true
      Boolean use_ensembl_annotation = false
      String gnomad_pop = "AF"
      Int mem_mb = 8000
      Int disk_space = 50
      Int cpu = 1
      String docker
      AnnotationRuntime runtime_params
    }

    String output_prefix = sample_id + ".chip_annotations"

    String ensembl_param = if (use_ensembl_annotation) then "--ensembl" else ""
    String gnomad_param = if (use_gnomad_genome) then "" else "--exome"
    String missing_af_param = if (treat_missing_as_rare) then "" else "--discard_missing_gnomad"
    String af_param = if (gnomad_pop == "AF") then "" else "--gnomad_population" + gnomad_pop

    command <<<
      annotate_chip \
        --sample ~{sample_id} \
        --vcf_header ~{input_vcf_header} \
        --chip_definitions ~{chip_mutations_csv} \
        --seq ~{seq_context_tsv} \
        --annovar ~{annovar_txt_input} \
        --annovar_function ~{annovar_variant_function} \
        --annovar_exonic_function ~{annovar_exonic_function} \
        --somaticism_transcripts ~{somaticism_filter_transcripts} \
        ~{ensembl_param} \
        ~{gnomad_param} \
        ~{missing_af_param} \
        ~{af_param} \
        --output_prefix ~{output_prefix}

      bcftools sort ~{output_prefix}.vcf > ~{output_prefix}.sorted.vcf
      bgzip -c ~{output_prefix}.sorted.vcf > ~{output_prefix}.sorted.vcf.gz
      tabix -s 1 -b 2 -e 2 ~{output_prefix}.sorted.vcf.gz

      bcftoolss sort ~{output_prefix}.merged.vcf > ~{output_prefix}.merged.sorted.vcf
      bgzip -c ~{output_prefix}.merged.sorted.vcf > ~{output_prefix}.merged.sorted.vcf.gz
      tabix -s 1 -b 2 -e 2 ~{output_prefix}.merged.sorted.vcf.gz
    >>>

    runtime {
      docker: docker
      bootDiskSizeGb: runtime_params.boot_disk_size
      memory: mem_mb + " MB"
      disks: "local-disk " + disk_space + " HDD"
      preemptible: runtime_params.preemptible
      maxRetries: runtime_params.max_retries
      cpu: cpu
    }

    output {
      File chip_vcf = output_prefix + ".merged.sorted.vcf.gz"
      File chip_vcf_idx = output_prefix + ".merged.sorted.vcf.gz.tbi"
      File chip_split_vcf = output_prefix + ".sorted.vcf.gz"
      File chip_split_vcf_idx = output_prefix + ".sorted.vcf.gz.tbi"
      File chip_csv = output_prefix + ".csv"
      File chip_rdata = output_prefix + ".RData"
    }
}

task FinaliseCHIPFilter {
  input {
      File input_vcf
      File input_vcf_idx
      File non_chip_gene_annotations_vcf
      File non_chip_gene_annotations_vf_idx
      File chip_gene_annotations_vcf
      File chip_gene_annotations_vcf_idx
      String sample_id
      Int mem_mb = 8000
      Int disk_space = 50
      Int cpu = 1
      String docker
      AnnotationRuntime runtime_params
    }

    String tmp_vcf = sample_id + ".tmp.vcf"
    String output_vcf = sample_id + ".chip.vcf"

    command <<<
      bcftools annotate \
        -a ~{non_chip_gene_annotations_vcf} \
        -c "=FILTER" \
        -o ~{tmp_vcf} \
        ~{input_vcf}
      bgzip -c ~{tmp_vcf} > ~{tmp_vcf}.gz
      tabix -f -s 1 -b 2 -e 2 ~{tmp_vcf}.gz
      bcftools annotate \
        -a ~{chip_gene_annotations_vcf} \
        -c "=FILTER,+INFO" \
        -o ~{output_vcf} \
        ~{tmp_vcf}.gz
      bgzip -c ${output_vcf} > ${output_vcf}.gz
      tabix -f -s 1 -b 2 -e 2 ${output_vcf}.gz
    >>>

    runtime {
      docker: docker
      bootDiskSizeGb: runtime_params.boot_disk_size
      memory: mem_mb + " MB"
      disks: "local-disk " + disk_space + " HDD"
      preemptible: runtime_params.preemptible
      maxRetries: runtime_params.max_retries
      cpu: cpu
    }

    output {
      File chip_vcf = output_vcf + ".gz"
      File chip_vcf_idx = output_vcf + ".gz.tbi"
    }
}
