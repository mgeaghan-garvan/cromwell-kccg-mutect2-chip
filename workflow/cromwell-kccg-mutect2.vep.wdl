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
#                                                                 #
# This pipeline has been developed for use by the Kinghorn        #
# Centre for Clinical Genomics and the Garvan Institute for       # 
# Medical Research.                                               #
#                                                                 #
# Author: Michael Geaghan (micgea)                                #
# Created: 2021/08/13                                             #
# =============================================================== #

import "cromwell-kccg-mutect2.wdl" as m2

workflow Mutect2CHIP_VEP {
    input {
        # input vcf
        File input_vcf
        File input_vcf_idx
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
        Int vep_tmp_disk = 100
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

    Int n_vep_cpus = if loftee then 1 else vep_cpu

    call m2.VEP {
        input:
            input_vcf = input_vcf,
            input_vcf_idx = input_vcf_idx,
            species = vep_species,
            assembly = vep_assembly,
            cache_archive = vep_cache_archive,
            fasta = ref_fasta,
            vep_docker = vep_docker,
            loftee_docker = loftee_docker,
            loftee = loftee,
            loftee_ancestor_fa = vep_loftee_ancestor_fa,
            loftee_ancestor_fai = vep_loftee_ancestor_fai,
            loftee_ancestor_gzi = vep_loftee_ancestor_gzi,
            loftee_conservation_sql = vep_loftee_conservation_sql,
            mem_mb = vep_mem,
            cpus = n_vep_cpus,
            vep_tmp_disk = vep_tmp_disk,
            runtime_params = standard_runtime
        }

    output {
        File? out_vep_vcf = VEP.output_vcf
    }
}