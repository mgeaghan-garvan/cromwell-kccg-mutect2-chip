version 1.0

# =============================================================== #
# extract_pass.wdl                                                #
#                                                                 #
# Extract PASS variants from a VCF with bcftools                  #
#                                                                 #
# This pipeline has been developed for use by the Kinghorn        #
# Centre for Clinical Genomics and the Garvan Institute for       # 
# Medical Research.                                               #
#                                                                 #
# Author: Michael Geaghan (micgea)                                #
# Created: 2024/10/28                                             #
# =============================================================== #

workflow ExtractPass {
    input {
        # input vcf
        File input_vcf
        File? input_vcf_idx
        String docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/bcftools:1.20"
        # runtime parameters
        Int disk_space = 100
        Int mem_mb = 4000
        Int preemptible = 2
        Int max_retries = 2
        Int boot_disk_size = 12
    }

    String sample_id = basename(basename(input_vcf, ".gz"), ".vcf")

    call ExtractPassVariants {
        input:
            input_vcf = input_vcf,
            input_vcf_idx = input_vcf_idx,
            sample_id = sample_id,
            docker = docker,
            disk_space = disk_space,
            mem_mb = mem_mb,
            preemptible = preemptible,
            max_retries = max_retries,
            boot_disk_size = boot_disk_size
    }

    output {
        File out_pass_variants = ExtractPassVariants.pass_variants
        File out_pass_variants_idx = ExtractPassVariants.pass_variants_idx
    }
}

task ExtractPassVariants {
    input {
      File input_vcf
      File? input_vcf_idx
      String sample_id
      String docker
      Int disk_space = 100
      Int mem_mb = 4000
      Int preemptible = 2
      Int max_retries = 2
      Int boot_disk_size = 12
    }

    String tabix = if (defined(input_vcf_idx)) then "FALSE" else "TRUE"

    command <<<
      set -euo pipefail

      # BGZIP and/or tabix if necessary
      BGZIP="$(basename ~{input_vcf} | sed -E -e 's|^.*\.([^\.]+)$|\1|g')"
      if [ ! "$BGZIP" == "gz" ]
      then
        bgzip -c ~{input_vcf} > ~{input_vcf}.gz
        tabix -s1 -b2 -e2 ~{input_vcf}.gz
        VCF=~{input_vcf}.gz
      elif [ "~{tabix}" == "TRUE" ]
      then
        tabix -s1 -b2 -e2 ~{input_vcf}
        VCF=~{input_vcf}
      fi

      # Extract PASS variants
      bcftools view -i 'FILTER=="PASS"' $VCF > ~{sample_id}.pass.vcf
      bgzip -c ~{sample_id}.pass.vcf > ~{sample_id}.pass.vcf.gz
      tabix -s1 -b2 -e2 ~{sample_id}.pass.vcf.gz
    >>>

    runtime {
      docker: docker
      bootDiskSizeGb: boot_disk_size
      memory: mem_mb + " MB"
      disks: "local-disk " + disk_space + " HDD"
      preemptible: preemptible
      maxRetries: max_retries
      cpu: 1
    }

    output {
      File pass_variants = sample_id + ".pass.vcf.gz"
      File pass_variants_idx = sample_id + ".pass.vcf.gz.tbi"
    }
}
