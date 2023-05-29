version 1.0

# =============================================================== #
# vep.wdl                                                         #
#                                                                 #
# This workflow annotates a VCF using VEP, and optionally with    #
# the LOFTEE VEP plugin.                                          #
#                                                                 #
# This pipeline has been developed for use by the Kinghorn        #
# Centre for Clinical Genomics and the Garvan Institute for       # 
# Medical Research.                                               #
#                                                                 #
# Author: Michael Geaghan (micgea)                                #
# Created: 2023/04/04                                             #
# =============================================================== #

struct VEPRuntime {
    Int max_retries
    Int preemptible
    Int disk
    Int boot_disk_size
}

workflow VEP {
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
        Boolean loftee = true
        File? vep_loftee_ancestor_fa
        File? vep_loftee_ancestor_fai
        File? vep_loftee_ancestor_gzi
        File? vep_loftee_conservation_sql
        String loftee_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/vep-loftee@sha256:c95b78bacef4c8d3770642138e6f28998a5034cfad3fbef5451d2303c8c795d3"  # :vep_103.1_loftee_1.0.3
        Int vep_mem = 32000
        Int vep_disk = 100
        Int vep_tmp_disk = 100
        Int vep_cpu = 1
        # runtime parameters
        Int? preemptible
        Int? max_retries
        Int boot_disk_size = 12
        # Use as a last resort to increase the disk given to every task in case of ill behaving data
        Int? emergency_extra_disk
    }

    Int preemptible_or_default = select_first([preemptible, 2])
    Int max_retries_or_default = select_first([max_retries, 2])

    # This is added to every task as padding, should increase if systematically you need more disk for every call
    Int disk_pad = 10 + select_first([emergency_extra_disk,0])

    VEPRuntime standard_runtime = {
        "max_retries": max_retries_or_default,
        "preemptible": preemptible_or_default,
        "disk": vep_disk + disk_pad,
        "boot_disk_size": boot_disk_size
    }

    Int n_vep_cpus = if loftee then 1 else vep_cpu

    call VEP_task {
        input:
            input_vcf = input_vcf,
            input_vcf_idx = input_vcf_idx,
            species = vep_species,
            assembly = vep_assembly,
            cache_archive = vep_cache_archive,
            fasta = ref_fasta,
            vep_docker = vep_docker,
            loftee = loftee,
            loftee_ancestor_fa = vep_loftee_ancestor_fa,
            loftee_ancestor_fai = vep_loftee_ancestor_fai,
            loftee_ancestor_gzi = vep_loftee_ancestor_gzi,
            loftee_conservation_sql = vep_loftee_conservation_sql,
            loftee_docker = loftee_docker,
            mem_mb = vep_mem,
            cpus = n_vep_cpus,
            vep_tmp_disk = vep_tmp_disk,
            runtime_params = standard_runtime
        }

    output {
        File? out_vep_vcf = VEP_task.output_vcf
    }
}

task VEP_task {
    input {
        # Need to be optional since the input file depends on whether realignment artifacts have been filtered
        File input_vcf
        File input_vcf_idx
        String species = "homo_sapiens"
        String assembly = "GRCh38"
        File cache_archive
        File fasta
        String vep_docker
        Boolean loftee = true
        File? loftee_ancestor_fa
        File? loftee_ancestor_fai
        File? loftee_ancestor_gzi
        File? loftee_conservation_sql
        String loftee_docker
        Boolean offline = true
        Int cpus = 1
        Int mem_mb = 32000
        Int? buffer_size
        Int loftee_buffer_size = 1
        Int vep_tmp_disk = 100
        VEPRuntime runtime_params
    }

    String output_options = "--vcf --no_stats"
    String input_basename = basename(basename(input_vcf, ".gz"), ".vcf")
    String output_suffix = ".vep.vcf"
    String output_file = input_basename + output_suffix
    String stats_file = input_basename + ".vep.html"
    String offline_options = if offline then "--offline" else ""

    # DNAnexus compatability: get the filename of all optional index files
    String loftee_ancestor_fai_def = if defined(loftee_ancestor_fai) then "defined" else "undefined"
    String loftee_ancestor_gzi_def = if defined(loftee_ancestor_gzi) then "defined" else "undefined"

    String tmp_dir = "/tmp"

    command {
        # DNAnexus compatability: echo optional index filenames to ensure they get localized
        OPT_VAR_DEFINED="~{loftee_ancestor_fai_def}"
        OPT_VAR_DEFINED="~{loftee_ancestor_gzi_def}"

        mkdir -p ~{tmp_dir}/.vep/cache
        tar -xzvf ~{cache_archive} -C ~{tmp_dir}/.vep/cache/
        # this seems necessary on GCP - running into permissions errors.
        # TODO: find a better solution
        cp ~{fasta} ~{tmp_dir}/ref_fasta.fasta
        vep \
            --dir_cache ~{tmp_dir}/.vep/cache \
            --dir_plugins /plugins/loftee-1.0.3 \
            -i ~{input_vcf} \
            --species ~{species} \
            --assembly ~{assembly} \
            ~{output_options} \
            -o ~{output_file} \
            --stats_file ~{stats_file} \
            ~{offline_options} \
            --cache \
            --fasta ~{tmp_dir}/ref_fasta.fasta \
            ~{if loftee then "" else "--fork " + cpus} \
            ~{if loftee then "--buffer_size " + loftee_buffer_size else if defined(buffer_size) then "--buffer_size " + select_first([buffer_size, 1]) else ""} \
            --no_progress \
            --everything \
            --pubmed \
            --hgvsg \
            --shift_hgvs 1 \
            ~{if loftee then "--plugin LoF,loftee_path:/plugins/loftee-1.0.3,human_ancestor_fa:" + loftee_ancestor_fa + ",conservation_file:" + loftee_conservation_sql else ""}
    }

    String docker_to_use = if (loftee && loftee_docker != "") then loftee_docker else vep_docker
    runtime {
        docker: docker_to_use
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: mem_mb + " MB"
        disks: ["local-disk ~{runtime_params.disk} HDD", "/tmp ~{vep_tmp_disk} HDD"]
        tmp_disk: runtime_params.disk
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: cpus
    }

    output {
        File? output_vcf = "~{input_basename}.vep.vcf"
    }
}