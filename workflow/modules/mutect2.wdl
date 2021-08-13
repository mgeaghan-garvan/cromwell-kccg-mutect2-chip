version 1.0

task M2 {
    input {
      File? intervals
      File ref_fasta
      File ref_fai
      File ref_dict
      File tumor_bam
      File tumor_bai
      File? normal_bam
      File? normal_bai
      File? pon
      File? pon_idx
      File? gnomad
      File? gnomad_idx
      String? m2_extra_args
      Boolean? make_bamout
      Boolean? run_ob_filter
      Boolean compress
      File? gga_vcf
      File? gga_vcf_idx
      File? variants_for_contamination
      File? variants_for_contamination_idx

      File? gatk_override

      # runtime
      String gatk_docker
      Int? mem
      Int? preemptible
      Int? max_retries
      Int? disk_space
      Int? cpu
      Boolean use_ssd = false
    }

    String output_vcf = "output" + if compress then ".vcf.gz" else ".vcf"
    String output_vcf_idx = output_vcf + if compress then ".tbi" else ".idx"

    String output_stats = output_vcf + ".stats"

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1000 else 3500
    Int command_mem = machine_mem - 500

    parameter_meta{
      intervals: {localization_optional: true}
      ref_fasta: {localization_optional: true}
      ref_fai: {localization_optional: true}
      ref_dict: {localization_optional: true}
      tumor_bam: {localization_optional: true}
      tumor_bai: {localization_optional: true}
      normal_bam: {localization_optional: true}
      normal_bai: {localization_optional: true}
      pon: {localization_optional: true}
      pon_idx: {localization_optional: true}
      gnomad: {localization_optional: true}
      gnomad_idx: {localization_optional: true}
      gga_vcf: {localization_optional: true}
      gga_vcf_idx: {localization_optional: true}
      variants_for_contamination: {localization_optional: true}
      variants_for_contamination_idx: {localization_optional: true}
    }

    command <<<
        set -e

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

        # We need to create these files regardless, even if they stay empty
        touch bamout.bam
        touch f1r2.tar.gz
        echo "" > normal_name.txt

        gatk --java-options "-Xmx~{command_mem}m" GetSampleName -R ~{ref_fasta} -I ~{tumor_bam} -O tumor_name.txt -encode
        tumor_command_line="-I ~{tumor_bam} -tumor `cat tumor_name.txt`"

        if [[ ! -z "~{normal_bam}" ]]; then
            gatk --java-options "-Xmx~{command_mem}m" GetSampleName -R ~{ref_fasta} -I ~{normal_bam} -O normal_name.txt -encode
            normal_command_line="-I ~{normal_bam} -normal `cat normal_name.txt`"
        fi

        gatk --java-options "-Xmx~{command_mem}m" Mutect2 \
            -R ~{ref_fasta} \
            $tumor_command_line \
            $normal_command_line \
            ~{"--germline-resource " + gnomad} \
            ~{"-pon " + pon} \
            ~{"-L " + intervals} \
            ~{"--alleles " + gga_vcf} \
            -O "~{output_vcf}" \
            ~{true='--bam-output bamout.bam' false='' make_bamout} \
            ~{true='--f1r2-tar-gz f1r2.tar.gz' false='' run_ob_filter} \
            ~{m2_extra_args}

        m2_exit_code=$?

        ### GetPileupSummaries

        # If the variants for contamination and the intervals for this scatter don't intersect, GetPileupSummaries
        # throws an error.  However, there is nothing wrong with an empty intersection for our purposes; it simply doesn't
        # contribute to the merged pileup summaries that we create downstream.  We implement this by with array outputs.
        # If the tool errors, no table is created and the glob yields an empty array.
        set +e

        if [[ ! -z "~{variants_for_contamination}" ]]; then
            gatk --java-options "-Xmx~{command_mem}m" GetPileupSummaries -R ~{ref_fasta} -I ~{tumor_bam} ~{"--interval-set-rule INTERSECTION -L " + intervals} \
                -V ~{variants_for_contamination} -L ~{variants_for_contamination} -O tumor-pileups.table

            if [[ ! -z "~{normal_bam}" ]]; then
                gatk --java-options "-Xmx~{command_mem}m" GetPileupSummaries -R ~{ref_fasta} -I ~{normal_bam} ~{"--interval-set-rule INTERSECTION -L " + intervals} \
                    -V ~{variants_for_contamination} -L ~{variants_for_contamination} -O normal-pileups.table
            fi
        fi

        # the script only fails if Mutect2 itself fails
        exit $m2_exit_code
    >>>

    runtime {
        docker: gatk_docker
        bootDiskSizeGb: 12
        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, 100]) + if use_ssd then " SSD" else " HDD"
        preemptible: select_first([preemptible, 10])
        maxRetries: select_first([max_retries, 0])
        cpu: select_first([cpu, 1])
    }

    output {
        File unfiltered_vcf = "~{output_vcf}"
        File unfiltered_vcf_idx = "~{output_vcf_idx}"
        File output_bamOut = "bamout.bam"
        String tumor_sample = read_string("tumor_name.txt")
        String normal_sample = read_string("normal_name.txt")
        File stats = "~{output_stats}"
        File f1r2_counts = "f1r2.tar.gz"
        Array[File] tumor_pileups = glob("*tumor-pileups.table")
        Array[File] normal_pileups = glob("*normal-pileups.table")
    }
}