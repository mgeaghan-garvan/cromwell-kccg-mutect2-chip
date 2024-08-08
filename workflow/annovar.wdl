version 1.0

# =============================================================== #
# annovar.wdl                                                     #
#                                                                 #
# This workflow annotates a VCF using Annovar.                    #
#                                                                 #
# This pipeline has been developed for use by the Kinghorn        #
# Centre for Clinical Genomics and the Garvan Institute for       # 
# Medical Research.                                               #
#                                                                 #
# Author: Michael Geaghan (micgea)                                #
# Created: 2023/04/04                                             #
# =============================================================== #

struct AnnotationRuntime {
    Int max_retries
    Int preemptible
    Int boot_disk_size
}

workflow Annovar {
    input {
        # input vcf
        File input_vcf
        # annovar settings
        Int annovar_mem_mb = 4000
        Int annovar_disk = 100
        Int annovar_tmp_disk = 200
        Int annovar_cpu = 1
        String annovar_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/annovar@sha256:842e9f88dd39999ee2129aeb992e8eced10ac2a33642d4b34d0f0c0254aa5035"  # :5.34.0
        File annovar_db_archive
        String ref_name = "hg38"
        String annovar_protocols = "refGene"
        String annovar_operations = "g"
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
    call Annovar_task {
        input:
            vcf_input = input_vcf,
            sample_id = sample_id,
            mem_mb = annovar_mem_mb,
            annovar_disk_space = annovar_disk,
            annovar_tmp_disk_space = annovar_tmp_disk,
            cpu = annovar_cpu,
            annovar_docker = annovar_docker,
            annovar_db_archive = annovar_db_archive,
            ref_name = ref_name,
            annovar_protocols = annovar_protocols,
            annovar_operations = annovar_operations,
            runtime_params = standard_runtime
    }

    output {
        File out_annovar_vcf = Annovar_task.annovar_output_file_vcf
        File out_annovar_table = Annovar_task.annovar_output_file_table
    }
}

task Annovar_task {
    input {
      File vcf_input
      String sample_id

      Int mem_mb = 4000
      Int annovar_disk_space = 100
      Int annovar_tmp_disk_space = 200
      Int cpu = 1
      String annovar_docker
      File annovar_db_archive

      String label = "annovar_out"

      String ref_name = "hg38"
      String annovar_protocols = "refGene"
      String annovar_operations = "g"
      String annovar_additional_arguments = ""

      AnnotationRuntime runtime_params
    }

    String file_prefix = sample_id + "." + label

    String tmp_dir = "/tmp"

    command <<<
      set -euo pipefail

      # Extract and flatten the annovar db archive
      mkdir ~{tmp_dir}/extract
      mkdir ~{tmp_dir}/annovar_db

      tar -xzvf ~{annovar_db_archive} -C ~{tmp_dir}/extract

      cd ~{tmp_dir}/annovar_db
      find ~{tmp_dir}/extract -type f | while read FILE; do ln -s $FILE; done
      cd -

      table_annovar.pl \
        ~{vcf_input} \
        ~{tmp_dir}/annovar_db \
        -buildver ~{default="hg38" ref_name} \
        -out ~{file_prefix} \
        -protocol ~{annovar_protocols} \
        -operation ~{annovar_operations} \
        -nastring . \
        -vcfinput \
        -polish \
        ~{annovar_additional_arguments}

      mv ~{file_prefix}.hg38_multianno.vcf ~{file_prefix}.hg38_multianno.vcf.bad
      gawk -v FS="\t" -v OFS="\t" '$0 ~ /^#/ { print $0 } $0 !~ /^#/ { $8 = gensub("^\.;", "", "g", $8); print $0 }' ~{file_prefix}.hg38_multianno.vcf.bad > ~{file_prefix}.hg38_multianno.vcf
    >>>

    runtime {
      docker: annovar_docker
      bootDiskSizeGb: runtime_params.boot_disk_size
      memory: mem_mb + " MB"
      disks: ["local-disk ~{annovar_disk_space} HDD", "/tmp ~{annovar_tmp_disk_space} HDD"]
      preemptible: runtime_params.preemptible
      maxRetries: runtime_params.max_retries
      cpu: cpu
    }

    output {
      File annovar_output_file_vcf = file_prefix + ".hg38_multianno.vcf"
      File annovar_output_file_table = file_prefix + ".hg38_multianno.txt"
      File? annovar_output_refgene_variant_function = file_prefix + ".refGene.variant_function"
      File? annovar_output_ensgene_variant_function = file_prefix + ".ensGene.variant_function"
      File? annovar_output_refgene_variant_exonic_function = file_prefix + ".refGene.exonic_variant_function"
      File? annovar_output_ensgene_variant_exonic_function = file_prefix + ".ensGene.exonic_variant_function"
    }
}
