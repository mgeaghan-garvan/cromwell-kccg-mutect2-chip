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

import "m2.wdl" as Mutect2
import "vep.wdl" as VEP
import "annovar.wdl" as Annovar
import "spliceai.wdl" as SpliceAI
import "chip.wdl" as CHIP

workflow Mutect2CHIP {
    input {
        # Mutect2 inputs
        File? intervals
        File ref_fasta
        File ref_fai
        File ref_dict
        File tumor_reads
        File tumor_reads_index
        File? normal_reads
        File? normal_reads_index
        File? pon
        File? pon_idx
        Int scatter_count = 10
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
        Boolean vep = false
        File? vep_cache_archive
        String vep_species = "homo_sapiens"
        String vep_assembly = "GRCh38"
        String vep_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/vep@sha256:bc6a74bf271adb1484ea769660c7b69f5eea033d3ba2e2947988e6c5f034f221"  # :release_103.1
        Boolean loftee = false
        File? vep_loftee_ancestor_fa
        File? vep_loftee_ancestor_fai
        File? vep_loftee_ancestor_gzi
        File? vep_loftee_conservation_sql
        String loftee_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/vep-loftee@sha256:c95b78bacef4c8d3770642138e6f28998a5034cfad3fbef5451d2303c8c795d3"  # :vep_103.1_loftee_1.0.3

        # Annovar settings
        Boolean annovar = true
        File? annovar_archive
        String annovar_assembly = "hg38"
        String annovar_protocols = "refGene"
        String annovar_operations = "g"
        String annovar_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/perl@sha256:842e9f88dd39999ee2129aeb992e8eced10ac2a33642d4b34d0f0c0254aa5035"  # :5.34.0

        # SpliceAI settings
        Boolean spliceai = false
        File? spliceai_annotation_file
        String spliceai_annotation_string = 'grch38'
        Int spliceai_max_dist = 50
        Boolean spliceai_mask = false
        String spliceai_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/spliceai@sha256:617682496a3f475c69ccdfe593156b79dd1ba21e02481ed1d0d8b740f3422530"  # :v1.3.1

        # Whitelist Filter settings
        Boolean run_chip_detection = true
        Boolean treat_missing_as_rare = true
        Boolean whitelist_genome = true
        Boolean whitelist_use_ensembl_annotation = false
        Boolean run_chip_on_unannotated_vcf = false
        String gnomad_pop = "AF"
        String whitelist_filter_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/whitelist_filter@sha256:3c4e9d0a650427a028814605423a39b02c621357043a5cff694fe671ee5a77ce"  # :latest

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
        Int m2_mem = 5000
        Int m2_cpu = 4
        Int learn_read_orientation_mem = 5000
        Int filter_alignment_artifacts_mem = 5000
        Int vep_mem = 32000
        Int vep_cpu = 1
        Int vep_tmp_disk = 100
        Int annovar_mem_mb = 4000
        Int annovar_disk = 100
        Int annovar_tmp_disk = 200
        Int spliceai_disk = 100
        Int spliceai_mem_mb = 16000
        Int spliceai_cpu = 4
        Int whitelist_mem_mb = 10000
        Int whitelist_disk = 300

        # Use as a last resort to increase the disk given to every task in case of ill behaving data
        Int? emergency_extra_disk

        # These are multipliers to multipler inputs by to make sure we have enough disk to accommodate for possible output sizes
        # Large is for Bams/WGS vcfs
        # Small is for metrics/other vcfs
        Float large_input_to_output_multiplier = 2.25
        Float small_input_to_output_multiplier = 2.0
    }

    # Run M2
    call Mutect2.Mutect2 as Mutect2_wf {
        input:
            intervals = intervals,
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            tumor_reads = tumor_reads,
            tumor_reads_index = tumor_reads_index,
            normal_reads = normal_reads,
            normal_reads_index = normal_reads_index,
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
            gatk_docker = gatk_docker,
            gatk_override = gatk_override,
            preemptible = preemptible,
            max_retries = max_retries,
            small_task_cpu = small_task_cpu,
            small_task_mem = small_task_mem,
            small_task_disk = small_task_disk,
            command_mem_padding = command_mem_padding,
            boot_disk_size = boot_disk_size,
            m2_mem = m2_mem,
            m2_cpu = m2_cpu,
            learn_read_orientation_mem = learn_read_orientation_mem,
            filter_alignment_artifacts_mem = filter_alignment_artifacts_mem,
            emergency_extra_disk = emergency_extra_disk,
            large_input_to_output_multiplier = large_input_to_output_multiplier,
            small_input_to_output_multiplier = small_input_to_output_multiplier
    }

    # Optionally run VEP
    if (vep && defined(vep_cache_archive)) {
        File vep_cache_archive_file = select_first([vep_cache_archive, "CACHE_FILE_NOT_SUPPLIED"])
        call VEP.VEP as VEP_wf {
            input:
                input_vcf = Mutect2_wf.filtered_vcf,
                input_vcf_idx = Mutect2_wf.filtered_vcf_idx,
                vep_species = vep_species,
                vep_assembly = vep_assembly,
                vep_cache_archive = vep_cache_archive_file,
                ref_fasta = ref_fasta,
                vep_docker = vep_docker,
                loftee = loftee,
                vep_loftee_ancestor_fa = vep_loftee_ancestor_fa,
                vep_loftee_ancestor_fai = vep_loftee_ancestor_fai,
                vep_loftee_ancestor_gzi = vep_loftee_ancestor_gzi,
                vep_loftee_conservation_sql = vep_loftee_conservation_sql,
                loftee_docker = loftee_docker,
                vep_mem = vep_mem,
                vep_disk = small_task_disk,
                vep_tmp_disk = vep_tmp_disk,
                vep_cpu = vep_cpu,
                preemptible = preemptible,
                max_retries = max_retries,
                boot_disk_size = boot_disk_size,
                emergency_extra_disk = emergency_extra_disk
        }
    }

    File annovar_input_vcf = select_first([VEP_wf.out_vep_vcf, Mutect2_wf.filtered_vcf])
    File annovar_archive_file = select_first([annovar_archive, "ANNOVAR_ARCHIVE_NOT_SUPPLIED"])

    # Optionally run Annovar
    if (annovar && defined(annovar_archive)) {
        call Annovar.Annovar as Annovar_wf {
            input:
                input_vcf = annovar_input_vcf,
                annovar_mem_mb = annovar_mem_mb,
                annovar_disk = annovar_disk,
                annovar_tmp_disk = annovar_tmp_disk,
                annovar_cpu = 1,
                annovar_docker = annovar_docker,
                annovar_archive = annovar_archive_file,
                ref_name = annovar_assembly,
                annovar_protocols = annovar_protocols,
                annovar_operations = annovar_operations,
                preemptible = preemptible,
                max_retries = max_retries,
                boot_disk_size = boot_disk_size
        }
    }

    File splieai_input_vcf = select_first([Annovar_wf.out_annovar_vcf, annovar_input_vcf])

    # Optionally run SpliceAI
    if (spliceai) {
        call SpliceAI.SpliceAI as SpliceAI_wf {
            input:
                input_vcf = splieai_input_vcf,
                ref_fasta = ref_fasta,
                spliceai_annotation_file = spliceai_annotation_file,
                spliceai_annotation_string = spliceai_annotation_string,
                spliceai_max_dist = spliceai_max_dist,
                spliceai_mask = spliceai_mask,
                spliceai_docker = spliceai_docker,
                spliceai_disk = spliceai_disk,
                spliceai_mem_mb = spliceai_mem_mb,
                spliceai_cpu = spliceai_cpu,
                preemptible = preemptible,
                max_retries = max_retries,
                boot_disk_size = boot_disk_size
        }
    }

    File chip_detection_input_vcf = if (run_chip_on_unannotated_vcf) then Mutect2_wf.filtered_vcf else select_first([SpliceAI_wf.spliceai_output_vcf, splieai_input_vcf])
    
    # Optionally run CHIP
    if (run_chip_detection && defined(annovar_archive)) {
        String tumor_sample_name = Mutect2_wf.tumor_sample
        call CHIP.CHIP as CHIP_wf {
            input:
                input_vcf = chip_detection_input_vcf,
                tumor_sample_name = tumor_sample_name,
                annovar_mem_mb = annovar_mem_mb,
                annovar_disk = annovar_disk,
                annovar_tmp_disk = annovar_tmp_disk,
                annovar_cpu = 1,
                annovar_docker = annovar_docker,
                annovar_archive = annovar_archive_file,
                whitelist_mem_mb = whitelist_mem_mb,
                whitelist_disk = whitelist_disk,
                whitelist_cpu = 1,
                treat_missing_as_rare = treat_missing_as_rare,
                whitelist_genome = whitelist_genome,
                whitelist_use_ensembl_annotation = whitelist_use_ensembl_annotation,
                gnomad_pop = gnomad_pop,
                whitelist_docker = whitelist_filter_docker,
                ref_fasta = ref_fasta,
                ref_name = annovar_assembly,
                preemptible = preemptible,
                max_retries = max_retries,
                boot_disk_size = boot_disk_size
        }
    }

    output {
        File filtered_vcf = Mutect2_wf.filtered_vcf
        File filtered_vcf_idx = Mutect2_wf.filtered_vcf_idx
        File filtering_stats = Mutect2_wf.filtering_stats
        File mutect_stats = Mutect2_wf.mutect_stats
        File? contamination_table = Mutect2_wf.contamination_table
        File? bamout = Mutect2_wf.bamout
        File? bamout_index = Mutect2_wf.bamout_index
        File? maf_segments = Mutect2_wf.maf_segments
        File? read_orientation_model_params = Mutect2_wf.read_orientation_model_params
        File? out_vep_vcf = VEP_wf.out_vep_vcf
        File? out_annovar_vcf = Annovar_wf.out_annovar_vcf
        File? out_annovar_table = Annovar_wf.out_annovar_table
        File? out_spliceai_vcf = SpliceAI_wf.spliceai_output_vcf
        File? out_whitelist_annovar_vcf = CHIP_wf.out_whitelist_annovar_vcf
        File? out_whitelist_annovar_table = CHIP_wf.out_whitelist_annovar_table
        File? out_whitelist_annovar_output_refgene_variant_function = CHIP_wf.out_whitelist_annovar_output_refgene_variant_function
        File? out_whitelist_annovar_output_ensgene_variant_function = CHIP_wf.out_whitelist_annovar_output_ensgene_variant_function
        File? out_whitelist_annovar_output_refgene_variant_exonic_function = CHIP_wf.out_whitelist_annovar_output_refgene_variant_exonic_function
        File? out_whitelist_annovar_output_ensgene_variant_exonic_function = CHIP_wf.out_whitelist_annovar_output_ensgene_variant_exonic_function
        File? out_whitelist_filter_output_csv = CHIP_wf.out_whitelist_filter_output_csv
        File? out_whitelist_filter_output_allvariants_csv = CHIP_wf.out_whitelist_filter_output_allvariants_csv
        File? out_whitelist_filter_output_allvariantsfiltered_csv = CHIP_wf.out_whitelist_filter_output_allvariantsfiltered_csv
        File? out_whitelist_filter_output_exonicsplicingvariants_csv = CHIP_wf.out_whitelist_filter_output_exonicsplicingvariants_csv
        File? out_whitelist_filter_output_chiptranscriptvariants_csv = CHIP_wf.out_whitelist_filter_output_chiptranscriptvariants_csv
        File? out_whitelist_filter_output_chiptranscriptvariantsfiltered_csv = CHIP_wf.out_whitelist_filter_output_chiptranscriptvariantsfiltered_csv
        File? out_whitelist_filter_output_putativefilter_csv = CHIP_wf.out_whitelist_filter_output_putativefilter_csv
    }
}
