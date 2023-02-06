version 1.0

# =============================================================== #
# cromwell-kccg-mutect2-chip                                      #
#                                                                 #
# This workflow detects mutations associated with CHIP            #
# (clonal haematopoiesis of indeterminate potential) using the    #
# GATK4 Mutect2 somatic variant calling pipeline and the          #
# Cromwell workflow engine.                                       #
#                                                                 #
# This script will call the cromwell pipeline on multiple         #
# samples given an input TSV file specifying the input BAM        #
# and BAI files.                                                  #
#                                                                 #
# This workflow has been adapted from the CHIP-detection-Mutect2  #
# workflow developed by Alex Bick. The Mutect2 somatic variant    #
# calling stage is adapted from the GATK best practices workflow, #
# and the CHIP detection stage is adapted from the                #
# Annovar Whitelist Filter developed by Charlie Condon.           #
#                                                                 #
# This pipeline has been developed for use by the Kinghorn        #
# Centre for Clinical Genomics and the Garvan Institute for       # 
# Medical Research.                                               #
#                                                                 #
# Author: Michael Geaghan (micgea)                                #
# Created: 2021/08/13                                             #
# =============================================================== #

#  Run Mutect 2 on a list of tumors or tumor-normal pairs
#
#  Description of inputs
#  intervals: genomic intervals
#  ref_fasta, ref_fai, ref_dict: reference genome, index, and dictionary
#  pon, pon_idx: optional panel of normals and index in vcf format containing known false positves
#  scatter_count: number of parallel jobs when scattering over intervals
#  gnomad, gnomad_idx: optional database of known germline variants, obtainable from http://gnomad.broadinstitute.org/downloads
#  variants_for_contamination, variants_for_contamination_idx: vcf of common variants with allele frequencies fo calculating contamination
#  run_orientation_bias_filter: if true, run the orientation bias filter post-processing step
#  pair_list: a tab-separated table with no header in the following format:
#   TUMOR_1_BAM</TAB>TUMOR_1_bai</TAB>NORMAL_1_BAM</TAB>NORMAL_1_bai
#   TUMOR_2_BAM</TAB>TUMOR_2_bai</TAB>NORMAL_2_BAM</TAB>NORMAL_2_bai
#   . . .
#  Tumor-only input is the same but without the columns for the normal:
#  TUMOR_1_BAM</TAB>TUMOR_1_bai
#  TUMOR_2_BAM</TAB>TUMOR_2_bai
#   . . .

import "cromwell-kccg-mutect2.wdl" as m2

workflow Mutect2CHIP_Multi {
    input {
        File? intervals
        File ref_fasta
        File ref_fai
        File ref_dict
    
        File pair_list
    
        File? pon
        File? pon_idx
        Int scatter_count
        File? gnomad
        File? gnomad_idx
        File? variants_for_contamination
        File? variants_for_contamination_idx
        File? realignment_index_bundle
        String? realignment_extra_args
        Boolean? run_orientation_bias_mixture_model_filter
        String? m2_extra_args
        String? m2_extra_filtering_args
        String? split_intervals_extra_args
        Boolean? make_bamout
        Boolean? compress_vcfs
        File? gga_vcf
        File? gga_vcf_idx
    
        # VEP settings
        String vep_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/vep@sha256:bc6a74bf271adb1484ea769660c7b69f5eea033d3ba2e2947988e6c5f034f221"  # :release_103.1
        String loftee_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/vep-loftee@sha256:c95b78bacef4c8d3770642138e6f28998a5034cfad3fbef5451d2303c8c795d3"  # :vep_103.1_loftee_1.0.3
        Boolean vep = true
        String vep_species = "homo_sapiens"
        String vep_assembly = "GRCh38"
        File? vep_cache_archive
        Boolean loftee = true
        File? vep_loftee_ancestor_fa
        File? vep_loftee_ancestor_fai
        File? vep_loftee_ancestor_gzi
        File? vep_loftee_conservation_sql
    
        # Annovar settings
        Boolean annovar = false
        File? annovar_archive
        String annovar_assembly = "hg38"
        String annovar_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/perl@sha256:1f35086e2ff48dace3b3edeaa2ad1faf1e44c0612e00f00ea0fc1830b576a261"  # :5.34.0
        String annovar_protocols = "cosmic70"
        String annovar_operations = "f"

        # Whitelist Filter settings
        Boolean run_chip_detection = true
        File? whitelist_filter_archive
        Boolean treat_missing_as_rare = true
        Boolean whitelist_genome = true
        Boolean whitelist_use_ensembl_annotation = false
        Boolean run_chip_on_unannotated_vcf = false
        String gnomad_pop = "AF"
        String whitelist_filter_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/whitelist_filter@sha256:1f1f83f8241f40fbd1f21b19e2ccbdc184984fd9ec0b0a7bdfa97b8a73fed8a4"  # :latest
    
        # Samtools settings
        String samtools_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/vep-loftee@sha256:c95b78bacef4c8d3770642138e6f28998a5034cfad3fbef5451d2303c8c795d3"  # same as loftee_docker
    
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
        Int c2b_mem = 6000
        Int m2_mem = 5000
        Int m2_cpu = 4
        Int learn_read_orientation_mem = 5000
        Int filter_alignment_artifacts_mem = 5000
        Int vep_mem = 32000
        Int vep_cpu = 4
        Int vep_tmp_disk = 100
        Int annovar_mem_mb = 4000
        Int annovar_disk = 100
        Int annovar_tmp_disk = 200
        Int whitelist_mem_mb = 10000
        Int whitelist_disk = 300
    
        # Use as a last resort to increase the disk given to every task in case of ill behaving data
        Int? emergency_extra_disk
    
        # These are multipliers to multipler inputs by to make sure we have enough disk to accommodate for possible output sizes
        # Large is for Bams/WGS vcfs
        # Small is for metrics/other vcfs
        Float large_input_to_output_multiplier = 2.25
        Float small_input_to_output_multiplier = 2.0
        Float cram_to_bam_multiplier = 6.0
    }
  
    Array[Array[String]] pairs = read_tsv(pair_list)

    scatter( row in pairs ) {
        #      If the condition is true, variables inside the 'if' block retain their values outside the block.
        #      Otherwise they are treated as null, which in WDL is equivalent to an empty optional
        if(length(row) == 4) {
            File normal_bam = row[2]
            File normal_bai = row[3]
        }

        call m2.Mutect2CHIP {
            input:
                intervals = intervals,
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                ref_dict = ref_dict,
                tumor_reads = row[0],
                tumor_reads_index = row[1],
                normal_reads = normal_bam,
                normal_reads_index = normal_bai,
                pon = pon,
                pon_idx = pon_idx,
                scatter_count = scatter_count,
                gnomad = gnomad,
                gnomad_idx = gnomad_idx,
                variants_for_contamination = variants_for_contamination,
                variants_for_contamination_idx = variants_for_contamination_idx,
                realignment_index_bundle = realignment_index_bundle,
                realignment_extra_args = realignment_extra_args,
                run_orientation_bias_mixture_model_filter = run_orientation_bias_mixture_model_filter,
                m2_extra_args = m2_extra_args,
                m2_extra_filtering_args = m2_extra_filtering_args,
                split_intervals_extra_args = split_intervals_extra_args,
                make_bamout = make_bamout,
                compress_vcfs = compress_vcfs,
                gga_vcf = gga_vcf,
                gga_vcf_idx = gga_vcf_idx,
                vep_docker = vep_docker,
                loftee_docker = loftee_docker,
                vep = vep,
                vep_species = vep_species,
                vep_assembly = vep_assembly,
                vep_cache_archive = vep_cache_archive,
                loftee = loftee,
                vep_loftee_ancestor_fa = vep_loftee_ancestor_fa,
                vep_loftee_ancestor_fai = vep_loftee_ancestor_fai,
                vep_loftee_ancestor_gzi = vep_loftee_ancestor_gzi,
                vep_loftee_conservation_sql = vep_loftee_conservation_sql,
                annovar = annovar,
                annovar_archive = annovar_archive,
                annovar_assembly = annovar_assembly,
                annovar_docker = annovar_docker,
                annovar_protocols = annovar_protocols,
                annovar_operations = annovar_operations,
                run_chip_detection = run_chip_detection,
                whitelist_filter_archive = whitelist_filter_archive,
                treat_missing_as_rare = treat_missing_as_rare,
                whitelist_genome = whitelist_genome,
                whitelist_use_ensembl_annotation = whitelist_use_ensembl_annotation,
                run_chip_on_unannotated_vcf = run_chip_on_unannotated_vcf,
                gnomad_pop = gnomad_pop,
                whitelist_filter_docker = whitelist_filter_docker,
                samtools_docker = samtools_docker,
                gatk_docker = gatk_docker,
                gatk_override = gatk_override,
                preemptible = preemptible,
                max_retries = max_retries,
                small_task_cpu = small_task_cpu,
                small_task_mem = small_task_mem,
                small_task_disk = small_task_disk,
                command_mem_padding = command_mem_padding,
                boot_disk_size = boot_disk_size,
                c2b_mem = c2b_mem,
                m2_mem = m2_mem,
                m2_cpu = m2_cpu,
                learn_read_orientation_mem = learn_read_orientation_mem,
                filter_alignment_artifacts_mem = filter_alignment_artifacts_mem,
                vep_mem = vep_mem,
                vep_cpu = vep_cpu,
                emergency_extra_disk = emergency_extra_disk,
                large_input_to_output_multiplier = large_input_to_output_multiplier,
                small_input_to_output_multiplier = small_input_to_output_multiplier,
                cram_to_bam_multiplier = cram_to_bam_multiplier,
                vep_tmp_disk = vep_tmp_disk,
                annovar_mem_mb = annovar_mem_mb,
                annovar_disk = annovar_disk,
                annovar_tmp_disk = annovar_tmp_disk,
                whitelist_mem_mb = whitelist_mem_mb,
                whitelist_disk = whitelist_disk
        }
    }

    output {
        Array[File] filtered_vcf = Mutect2CHIP.filtered_vcf 
        Array[File] filtered_vcf_idx = Mutect2CHIP.filtered_vcf_idx 
        Array[File] filtering_stats = Mutect2CHIP.filtering_stats 
        Array[File] mutect_stats = Mutect2CHIP.mutect_stats 
        Array[File?] contamination_table = Mutect2CHIP.contamination_table 
        Array[File?] bamout = Mutect2CHIP.bamout 
        Array[File?] bamout_index = Mutect2CHIP.bamout_index 
        Array[File?] maf_segments = Mutect2CHIP.maf_segments 
        Array[File?] read_orientation_model_params = Mutect2CHIP.read_orientation_model_params 
        Array[File?] out_vep_vcf = Mutect2CHIP.out_vep_vcf 
        Array[File?] out_annovar_vcf = Mutect2CHIP.out_annovar_vcf
        Array[File?] out_annovar_table = Mutect2CHIP.out_annovar_table
        Array[File?] out_whitelist_annovar_vcf = Mutect2CHIP.out_whitelist_annovar_vcf
        Array[File?] out_whitelist_annovar_table = Mutect2CHIP.out_whitelist_annovar_table
        Array[File?] out_whitelist_annovar_output_refgene_variant_function = Mutect2CHIP.out_whitelist_annovar_output_refgene_variant_function
        Array[File?] out_whitelist_annovar_output_ensgene_variant_function = Mutect2CHIP.out_whitelist_annovar_output_ensgene_variant_function
        Array[File?] out_whitelist_annovar_output_refgene_variant_exonic_function = Mutect2CHIP.out_whitelist_annovar_output_refgene_variant_exonic_function
        Array[File?] out_whitelist_annovar_output_ensgene_variant_exonic_function = Mutect2CHIP.out_whitelist_annovar_output_ensgene_variant_exonic_function
        Array[File?] out_whitelist_filter_output_allvariants_csv = Mutect2CHIP.out_whitelist_filter_output_allvariants_csv
        Array[File?] out_whitelist_filter_output_allvariantsfiltered_csv = Mutect2CHIP.out_whitelist_filter_output_allvariantsfiltered_csv
        Array[File?] out_whitelist_filter_output_exonicsplicingvariants_csv = Mutect2CHIP.out_whitelist_filter_output_exonicsplicingvariants_csv
        Array[File?] out_whitelist_filter_output_chiptranscriptvariants_csv = Mutect2CHIP.out_whitelist_filter_output_chiptranscriptvariants_csv
        Array[File?] out_whitelist_filter_output_chiptranscriptvariantsfiltered_csv = Mutect2CHIP.out_whitelist_filter_output_chiptranscriptvariantsfiltered_csv
        Array[File?] out_whitelist_filter_output_putativefilter_csv = Mutect2CHIP.out_whitelist_filter_output_putativefilter_csv
    }
}
