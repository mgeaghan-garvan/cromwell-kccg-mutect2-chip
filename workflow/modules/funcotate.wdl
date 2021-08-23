version 1.0

import "runtime.wdl" as RT

task Funcotate {
     input {
       File ref_fasta
       File ref_fai
       File ref_dict
       File input_vcf
       File input_vcf_idx
       String reference_version
       String output_file_base_name
       String output_format
       Boolean compress
       Boolean use_gnomad
       # This should be updated when a new version of the data sources is released
       # TODO: Make this dynamically chosen in the command.
       File? data_sources_tar_gz = "gs://broad-public-datasets/funcotator/funcotator_dataSources.v1.6.20190124s.tar.gz"
       String? control_id
       String? case_id
       String? sequencing_center
       String? sequence_source
       String? transcript_selection_mode
       File? transcript_selection_list
       Array[String]? annotation_defaults
       Array[String]? annotation_overrides
       Array[String]? funcotator_excluded_fields
       Boolean? filter_funcotations
       File? interval_list

       String? extra_args

       # ==============
       Runtime runtime_params
       Int? disk_space   #override to request more disk than default small task params

       # You may have to change the following two parameter values depending on the task requirements
       Int default_ram_mb = 3000
       # WARNING: In the workflow, you should calculate the disk space as an input to this task (disk_space_gb).  Please see [TODO: Link from Jose] for examples.
       Int default_disk_space_gb = 100
     }

     # ==============
     # Process input args:
     String output_maf = output_file_base_name + ".maf"
     String output_maf_index = output_maf + ".idx"
     String output_vcf = output_file_base_name + if compress then ".vcf.gz" else ".vcf"
     String output_vcf_idx = output_vcf +  if compress then ".tbi" else ".idx"
     String output_file = if output_format == "MAF" then output_maf else output_vcf
     String output_file_index = if output_format == "MAF" then output_maf_index else output_vcf_idx
     String transcript_selection_arg = if defined(transcript_selection_list) then " --transcript-list " else ""
     String annotation_def_arg = if defined(annotation_defaults) then " --annotation-default " else ""
     String annotation_over_arg = if defined(annotation_overrides) then " --annotation-override " else ""
     String filter_funcotations_args = if defined(filter_funcotations) && (filter_funcotations) then " --remove-filtered-variants " else ""
     String excluded_fields_args = if defined(funcotator_excluded_fields) then " --exclude-field " else ""
     String interval_list_arg = if defined(interval_list) then " -L " else ""
     String extra_args_arg = select_first([extra_args, ""])

     String dollar = "$"

     parameter_meta{
      ref_fasta: {localization_optional: true}
      ref_fai: {localization_optional: true}
      ref_dict: {localization_optional: true}
      input_vcf: {localization_optional: true}
      input_vcf_idx: {localization_optional: true}
     }

     command <<<
         set -e
         export GATK_LOCAL_JAR=~{default="/gatk/gatk.jar" runtime_params.gatk_override}

         # Extract our data sources:
         echo "Extracting data sources zip file..."
         mkdir datasources_dir
         tar zxvf ~{data_sources_tar_gz} -C datasources_dir --strip-components 1
         DATA_SOURCES_FOLDER="$PWD/datasources_dir"

         # Handle gnomAD:
         if ~{use_gnomad} ; then
             echo "Enabling gnomAD..."
             for potential_gnomad_gz in gnomAD_exome.tar.gz gnomAD_genome.tar.gz ; do
                 if [[ -f ~{dollar}{DATA_SOURCES_FOLDER}/~{dollar}{potential_gnomad_gz} ]] ; then
                     cd ~{dollar}{DATA_SOURCES_FOLDER}
                     tar -zvxf ~{dollar}{potential_gnomad_gz}
                     cd -
                 else
                     echo "ERROR: Cannot find gnomAD folder: ~{dollar}{potential_gnomad_gz}" 1>&2
                     false
                 fi
             done
         fi

         # Run Funcotator:
         gatk --java-options "-Xmx~{runtime_params.command_mem}m" Funcotator \
             --data-sources-path $DATA_SOURCES_FOLDER \
             --ref-version ~{reference_version} \
             --output-file-format ~{output_format} \
             -R ~{ref_fasta} \
             -V ~{input_vcf} \
             -O ~{output_file} \
             ~{interval_list_arg} ~{default="" interval_list} \
             --annotation-default normal_barcode:~{default="Unknown" control_id} \
             --annotation-default tumor_barcode:~{default="Unknown" case_id} \
             --annotation-default Center:~{default="Unknown" sequencing_center} \
             --annotation-default source:~{default="Unknown" sequence_source} \
             ~{"--transcript-selection-mode " + transcript_selection_mode} \
             ~{transcript_selection_arg}~{default="" sep=" --transcript-list " transcript_selection_list} \
             ~{annotation_def_arg}~{default="" sep=" --annotation-default " annotation_defaults} \
             ~{annotation_over_arg}~{default="" sep=" --annotation-override " annotation_overrides} \
             ~{excluded_fields_args}~{default="" sep=" --exclude-field " funcotator_excluded_fields} \
             ~{filter_funcotations_args} \
             ~{extra_args_arg}
         # Make sure we have a placeholder index for MAF files so this workflow doesn't fail:
         if [[ "~{output_format}" == "MAF" ]] ; then
            touch ~{output_maf_index}
         fi
     >>>

    runtime {
        docker: runtime_params.gatk_docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        mem_mb: runtime_params.machine_mem
        disks: "local-disk " + select_first([disk_space, runtime_params.disk]) + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

     output {
         File funcotated_output_file = "~{output_file}"
         File funcotated_output_file_index = "~{output_file_index}"
     }
}