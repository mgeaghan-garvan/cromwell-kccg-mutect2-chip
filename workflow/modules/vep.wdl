version 1.0

task VEP {
    input {
        File input_vcf
        File input_vcf_idx
        String species = "homo_sapiens"
        String assembly = "GRCh38"
        Boolean vcf_out = true
        File cache_dir = "${HOME}/.vep/"
        Boolean offline = true
        Int cpus = 4
        Int mem_mb = 6000
    }

    String output_options = if vcf_out then "--vcf --no_stats" else "--tab"
    String input_basename = basename(basename(input_vcf, ".gz"), ".vcf")
    String output_suffix = if vcf_out then ".vep.vcf" else ".vep.txt"
    String output_file = input_basename + output_suffix
    String stats_file = input_basename + ".vep.html"
    String offline_options = if offline then "--offline" else ""

    command {
        vep \
            --dir ~{cache_dir} \
            -i ~{input_vcf} \
            --species ~{species} \
            --assembly ~{assembly} \
            ~{output_options} \
            -o ~{output_file} \
            --stats_file ~{stats_file} \
            ~{offline_options} \
            --cache \
            --fork ~{cpus} \
            --no_progress \
            --everything \
            --pubmed \
            --hgvsg \
            --shift_hgvs 1
    }

    runtime {
        docker: "ensemblorg/ensembl-vep:release_103.1"
        mem_mb: mem_mb
        cpu: cpus
    }

    output {
        File? output_vcf = "~{input_basename}.vep.vcf"
        File? output_tab = "~{input_basename}.vep.txt"
    }
}