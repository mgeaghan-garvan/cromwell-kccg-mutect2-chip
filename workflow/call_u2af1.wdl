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


struct Runtime {
    String docker
    Int max_retries
    Int preemptible
    Int cpu
    Int mem_mb
    Int disk
    Int boot_disk_size
}

workflow CallU2AF1 {
    input {
        # Inputs
        File ref_fasta
        File ref_fai
        File tumor_reads
        File tumor_reads_index
        File u2af1_regions_file
        File mutect2_output_vcf
        File mutect2_output_vcf_index
        Boolean compress_vcfs = false
        
        # Runtime options
        String docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/chip_pre_post_filter:latest"
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
    Int mutect2_output_vcf_size = ceil(size(mutect2_output_vcf, "GB") + size(mutect2_output_vcf_index, "GB"))
    Int output_size = mutect2_output_vcf_size * 2
    Int total_size = (ref_size + tumor_reads_size + mutect2_output_vcf_size + output_size) * 1.5
    Int disk_size = if (total_size > disk) then total_size else disk

    Runtime standard_runtime = {
        "docker": docker,
        "max_retries": max_retries,
        "preemptible": preemptible,
        "cpu": cpu,
        "mem_mb": mem_mb,
        "disk": disk_size,
        "boot_disk_size": boot_disk_size
    }


    

    output {
        File u2af1_vcf = FilterU2AF1Vcf.filtered_vcf
        File u2af1_vcf_idx = FilterU2AF1Vcf.filtered_vcf_idx
        File merged_vcf = MergeU2AF1Vcf.merged_vcf
        File merged_vcf_idx = MergeU2AF1Vcf.merged_vcf_idx
    }
}

# ================ #
# TASK DEFINITIONS #
# ================ #

task U2AF1Pileup {
    input {
      File tumor_reads
      File tumor_reads_index
      File u2af1_regions_file
      File ref_fasta
      File ref_fai

      # runtime
      Runtime runtime_params
    }

    String sample_basename = basename(basename(tumor_reads, ".bam"),".cram")

    # The U2AF1 regions file is expected to be a TSV file with the following columns:
    # Chromosome, Start (1-based), End (1-based), Duplication start (1-based), Duplication end (1-based)

    command <<<
        # Generate a BED file from the U2AF1 regions file
        U2AF1_BED="$(basename ~{u2af1_regions_file}).bed"
        awk -v FS="\t" -v OFS="\t" '{
            print $1, $2 - 1, $3;
            print $1, $4 - 1, $5;
        }' ~{u2af1_regions_file} | bedtools sort > $U2AF1_BED

        # Run bcftools mpileup and call on the U2AF1 regions
        bcftools mpileup \
            -R $U2AF1_BED \
            -a "FORMAT/AD,FORMAT/DP" \
            -d 8000 \
            -f ~{ref_fasta} \
            ~{tumor_reads} | \
        bcftools call -mv -Ov -o ~{sample_basename}.u2af1.pileup.vcf
    >>>

    runtime {
        docker: runtime_params.docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.mem_mb + " MB"
        disks: "local-disk " + runtime_params.disk + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    output {
        File pileup_vcf = "~{sample_basename}.u2af1.pileup.vcf"
    }
}

task FilterU2AF1Vcf {
    input {
      File pileup_vcf
      File u2af1_regions_file

      # runtime
      Runtime runtime_params
    }

    command <<<
        R --vanilla <<EOF
            library(tidyverse)
            regions <- read_tsv("~{u2af1_regions_file}", col_names = c("chrom", "start", "end", "dup_start", "dup_end"))
            vcf <- read_tsv("~{pileup_vcf}", comment = "#", col_names = c("chrom", "pos", "id", "ref", "alt", "qual", "filter", "info", "format", "sample"), guess_max = Inf)
        EOF
    >>>

    runtime {
        docker: runtime_params.docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.mem_mb + " MB"
        disks: "local-disk " + runtime_params.disk + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    output {
        File filtered_vcf = ".u2af1.pileup.filtered.vcf.gz"
        File filtered_vcf_idx = ".u2af1.pileup.filtered.vcf.gz.tbi"
    }
}

task MergeU2AF1Vcf {
    input {
      File u2af1_vcf
      File u2af1_vcf_index
      File mutect2_vcf
      File mutect2_vcf_index

      # runtime
      Runtime runtime_params
    }

    command <<<

    >>>

    runtime {
        docker: runtime_params.docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.mem_mb + " MB"
        disks: "local-disk " + runtime_params.disk + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    output {
        File merged_vcf = ".u2af1.vcf.gz"
        File merged_vcf_idx = ".u2af1.vcf.gz.tbi"
    }
}
