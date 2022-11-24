version 1.0

# =============================================================== #
# cromwell-kccg-mutect2-chip                                      #
#                                                                 #
# This workflow detects mutations associated with CHIP            #
# (clonal haematopoiesis of indeterminate potential) using the    #
# GATK4 Mutect2 somatic variant calling pipeline and the          #
# Cromwell workflow engine.                                       #
#                                                                 #
# This script will call the VEP annotation pipeline on multiple   #
# samples given an input TSV file specifying the input VCF        #
# files.                                                          #
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

import "cromwell-kccg-mutect2.vep.wdl" as m2v

workflow Mutect2CHIP_VEP_Multi {
    input {
        # input vcf
        File input_vcf_list
        # VEP settings
        String vep_species = "homo_sapiens"
        String vep_assembly = "GRCh38"
        File vep_cache_archive
        File ref_fasta
        String vep_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/vep@sha256:bc6a74bf271adb1484ea769660c7b69f5eea033d3ba2e2947988e6c5f034f221"  # :release_103.1
        String loftee_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/vep-loftee@sha256:c95b78bacef4c8d3770642138e6f28998a5034cfad3fbef5451d2303c8c795d3"  # :vep_103.1_loftee_1.0.3
        Boolean loftee = true
        File? vep_loftee_ancestor_fa
        File? vep_loftee_ancestor_fai
        File? vep_loftee_ancestor_gzi
        File? vep_loftee_conservation_sql
        Int vep_mem = 32000
        Int vep_cpu = 1
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

    Array[Array[String]] input_vcfs = read_tsv(input_vcf_list)

    scatter( row in input_vcfs ) {
        File input_vcf = row[0]
        File input_vcf_idx = row[1]

        call m2v.Mutect2CHIP_VEP {
            input:
                input_vcf = input_vcf,
                input_vcf_idx = input_vcf_idx,
                vep_species = vep_species,
                vep_assembly = vep_assembly,
                vep_cache_archive = vep_cache_archive,
                ref_fasta = ref_fasta,
                vep_docker = vep_docker,
                loftee_docker = loftee_docker,
                loftee = loftee,
                vep_loftee_ancestor_fa = vep_loftee_ancestor_fa,
                vep_loftee_ancestor_fai = vep_loftee_ancestor_fai,
                vep_loftee_ancestor_gzi = vep_loftee_ancestor_gzi,
                vep_loftee_conservation_sql = vep_loftee_conservation_sql,
                vep_mem = vep_mem,
                vep_cpu = vep_cpu,
                preemptible = preemptible,
                max_retries = max_retries,
                small_task_cpu = small_task_cpu,
                small_task_mem = small_task_mem,
                small_task_disk = small_task_disk,
                command_mem_padding = command_mem_padding,
                boot_disk_size = boot_disk_size,
                emergency_extra_disk = emergency_extra_disk,
                use_tmp_dir = use_tmp_dir
        }
    }

    output {
        Array[File?] out_vep_vcf = Mutect2CHIP_VEP.out_vep_vcf
    }
}