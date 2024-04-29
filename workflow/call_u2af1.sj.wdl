version 1.0

# =================================================================================== #
# call_u2af1.wdl                                                                      #
#                                                                                     #
# This workflow is designed to call somatic CHIP mutations in the U2AF1 locus in the  #
# GRCh38 build of the human genome. This build contains a duplicated sequence of this #
# gene at chr21:6484623-6499248, which in Ensembl is annotated as the artifact        #
# "U2 small nuclear RNA auxiliary factor 1-like 5".                                   #
#                                                                                     #
# This script uses bcftools mpileup and bcftools call to generate a VCF of known CHIP #
# mutations in this region, then merge the calls into a single locus at the correct   #
# position of U2AF1 - chr21:43092956-43107570. It will then take the output VCF from  #
# Mutect2, mask out both loci, and merge in the new updated U2AF1 calls. It will also #
# run a simple filter to FAIL any multi-allelic mutations, similar to Mutect2.        #
#                                                                                     #
# It is recommended that downstream CHIP mutation calling classes these U2AF1         #
# mutations as putative, given the different calling method. The default chip         #
# mutation definition files in this repository have been updates accordingly.         #
#                                                                                     #
# This pipeline has been developed for use by the Kinghorn                            #
# Centre for Clinical Genomics and the Garvan Institute for                           # 
# Medical Research.                                                                   #
#                                                                                     #
# Author: Michael Geaghan (micgea)                                                    #
# Created: 2023/04/04                                                                 #
# =================================================================================== #

import "call_u2af1.wdl" as U2AF1

workflow CallU2AF1SingleJob {
    input {
        # Inputs
        File ref_fasta
        File ref_fai
        File tumor_reads
        File tumor_reads_index
        File u2af1_regions_file
        
        # Runtime options
        String pileup_docker
        Int preemptible = 2
        Int max_retries = 2
        Int cpu = 4
        Int mem_mb = 4000
        Int disk = 100
        Int boot_disk_size = 12
    }

    # Disk sizes used for dynamic sizing
    Int ref_size = ceil(size(ref_fasta, "GB") + size(ref_fai, "GB"))
    Int tumor_reads_size = ceil(size(tumor_reads, "GB") + size(tumor_reads_index, "GB"))
    Int total_size = ceil((ref_size + tumor_reads_size) * 1.5)
    Int disk_size = if (total_size > disk) then total_size else disk

    U2AF1Runtime standard_runtime = {
        "max_retries": max_retries,
        "preemptible": preemptible,
        "cpu": cpu,
        "mem_mb": mem_mb,
        "disk": disk_size,
        "boot_disk_size": boot_disk_size
    }

    call U2AF1.U2AF1Pileup {
        input:
            tumor_reads = tumor_reads,
            tumor_reads_index = tumor_reads_index,
            u2af1_regions_file = u2af1_regions_file,
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            docker = pileup_docker,
            runtime_params = standard_runtime
    }

    output {
        File u2af1_vcf = U2AF1Pileup.pileup_vcf
        File u2af1_tsv = U2AF1Pileup.pileup_tsv
        File u2af1_merged_tsv = U2AF1Pileup.pileup_merged_tsv
    }
}
