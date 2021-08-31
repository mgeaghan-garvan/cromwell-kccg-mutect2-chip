version 1.0

task VEP {
    input {
        # Need to be optional since the input file depends on whether realignment artifacts have been filtered
        File input_vcf
        File input_vcf_idx
        String species = "homo_sapiens"
        String assembly = "GRCh38"
        Boolean vcf_out = true
        String cache_dir
        Boolean loftee = true
        File? loftee_ancestor_fa
        File? loftee_ancestor_fai
        File? loftee_ancestor_gzi
        File? loftee_conservation_sql
        Boolean offline = true
        Int cpus = 1
        Int mem_mb = 32000
        Int? buffer_size
        Int loftee_buffer_size = 1
    }

    String output_options = if vcf_out then "--vcf --no_stats" else "--tab"
    String input_basename = basename(basename(input_vcf, ".gz"), ".vcf")
    String output_suffix = if vcf_out then ".vep.vcf" else ".vep.txt"
    String output_file = input_basename + output_suffix
    String stats_file = input_basename + ".vep.html"
    String offline_options = if offline then "--offline" else ""

    command {
        vep \
            --dir_cache /cache \
            --dir_plugins /plugins/loftee-1.0.3 \
            -i ~{input_vcf} \
            --species ~{species} \
            --assembly ~{assembly} \
            ~{output_options} \
            -o ~{output_file} \
            --stats_file ~{stats_file} \
            ~{offline_options} \
            --cache \
            ~{if loftee then "" else "--fork " + cpus} \
            ~{if loftee then "--buffer_size " + loftee_buffer_size else if defined(buffer_size) then "--buffer_size " + select_first([buffer_size, 1]) else ""} \
            --no_progress \
            --everything \
            --pubmed \
            --hgvsg \
            --shift_hgvs 1 \
            ~{if loftee then "--plugin LoF,loftee_path:/plugins/loftee-1.0.3,human_ancestor_fa:" + loftee_ancestor_fa + ",conservation_file:" + loftee_conservation_sql else ""}
    }

    String cache_dir_bind = "--bind ~{cache_dir}:/cache"

    runtime {
        docker: "ensemblorg/ensembl-vep:release_103.1"
        singularity_local: "/home/micgea/singularity/vep_loftee.sif"
        mem_mb: mem_mb
        cpu: cpus
        bind_cmd: cache_dir_bind
    }

    output {
        File? output_vcf = "~{input_basename}.vep.vcf"
        File? output_tab = "~{input_basename}.vep.txt"
    }
}