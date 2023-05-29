library(dplyr)

create_annot_table <- function(vars_df, vars_hf_hp_df, vars_ann_som_df, vars_chip_filtered_put_df, vars_annovar_vcf_ids_df, transcript, aachange, genedetail, func, exonic_func, refseq_ensembl_suffix) {    
  # Master table
  vcf_annot <- unique(vars_df[c(
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "FILTER"
  )])
  
  vcf_annot$INFO.start <- paste("CHIP_DATE", Sys.Date(), sep = "=")
  vcf_annot$INFO.end <- "CHIP_END"
  
  # ----- Initial filter table ----- #
  tmp_vcf_annot_vars_hf_hp <- vars_hf_hp_df[c(
    "Chr",
    "Start",
    "End",
    "Ref",
    "Alt",
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "FILTER",
    "VAF",
    "VAF_M2",
    "HARD_FILTER",
    "VAR_SEQ_CONTEXT",
    "VAR_ALT_SEQ_CONTEXT",
    "HOMOPOLYMER_FILTER"
  )]
  
  tmp_vcf_annot_vars_hf_hp$vcf_idx <- vars_annovar_vcf_ids_df[
    paste(
      tmp_vcf_annot_vars_hf_hp$Chr,
      tmp_vcf_annot_vars_hf_hp$Start,
      tmp_vcf_annot_vars_hf_hp$End,
      tmp_vcf_annot_vars_hf_hp$Ref,
      tmp_vcf_annot_vars_hf_hp$Alt,
      sep = ":"
    ),
    "idx"
  ]
  
  tmp_vcf_annot_vars_hf_hp <- tmp_vcf_annot_vars_hf_hp[c(
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "vcf_idx",
    "FILTER",
    "VAF",
    "VAF_M2",
    "HARD_FILTER",
    "VAR_SEQ_CONTEXT",
    "VAR_ALT_SEQ_CONTEXT",
    "HOMOPOLYMER_FILTER"
  )]
  
  tmp_vcf_annot_vars_hf_hp <- tmp_vcf_annot_vars_hf_hp[
    order(
      tmp_vcf_annot_vars_hf_hp$CHROM,
      tmp_vcf_annot_vars_hf_hp$POS,
      tmp_vcf_annot_vars_hf_hp$REF,
      tmp_vcf_annot_vars_hf_hp$ALT,
      tmp_vcf_annot_vars_hf_hp$vcf_idx
    ),
  ]
  
  tmp_vcf_annot_vars_hf_hp$INFO <- paste(
    paste("VAF", tmp_vcf_annot_vars_hf_hp$VAF, sep = "="),
    paste("VAF_M2", tmp_vcf_annot_vars_hf_hp$VAF_M2, sep = "="),
    paste("REF_SEQ_CTX", tmp_vcf_annot_vars_hf_hp$VAR_SEQ_CONTEXT, sep = "="),
    paste("ALT_SEQ_CTX", tmp_vcf_annot_vars_hf_hp$VAR_ALT_SEQ_CONTEXT, sep = "="),
    sep = ";"
  )
  
  tmp_vcf_annot_vars_hf_hp <- tmp_vcf_annot_vars_hf_hp %>%
    group_by(CHROM, POS, REF, ALT, FILTER) %>%
    mutate(
      hf = paste0(HARD_FILTER, collapse = ","),
      hp = paste0(HOMOPOLYMER_FILTER, collapse = ",")
    ) %>%
    summarise(
      FILTER.HARD = first(hf),
      FILTER.HOMOPOLYMER = first(hp),
      INFO.HFHP = first(INFO)
    ) %>%
    mutate(
      INFO.HARD_FILTER = paste("HARD_FILTER", FILTER.HARD, sep = "="),
      INFO.HOMOPOLYMER_FILTER = paste("HOMOPOLYMER_FILTER", FILTER.HOMOPOLYMER, sep = "=")
    )
  
  
  # ----- Annotation and somaticism filter table ----- #
  tmp_vcf_annot_vars_ann_som <- vars_ann_som_df[c(
    "Chr",
    "Start",
    "End",
    "Ref",
    "Alt",
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "FILTER",
    transcript,
    aachange,
    genedetail,
    func,
    exonic_func,
    "SOMATICISM_FILTER"
  )]
  
  colnames(tmp_vcf_annot_vars_ann_som) <- c(
    "Chr",
    "Start",
    "End",
    "Ref",
    "Alt",
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "FILTER",
    "Transcript",
    "AAChange",
    "GeneDetail",
    "Func",
    "ExonicFunc",
    "SOMATICISM_FILTER"
  )
  
  tmp_vcf_annot_vars_ann_som$vcf_idx <- vars_annovar_vcf_ids_df[
    paste(
      tmp_vcf_annot_vars_ann_som$Chr,
      tmp_vcf_annot_vars_ann_som$Start,
      tmp_vcf_annot_vars_ann_som$End,
      tmp_vcf_annot_vars_ann_som$Ref,
      tmp_vcf_annot_vars_ann_som$Alt,
      sep = ":"
    ),
    "idx"
  ]
  
  tmp_vcf_annot_vars_ann_som <- tmp_vcf_annot_vars_ann_som %>%
    group_by(Chr, Start, End, Ref, Alt) %>%
    mutate(
      tr = paste0(Transcript, collapse = "|"),
      aa = paste0(AAChange, collapse = "|"),
      gd = paste0(GeneDetail, collapse = "|"),
      fc = paste0(gsub("\\s+", "_", Func, perl = TRUE), collapse = "|"),
      ef = paste0(gsub("\\s+", "_", ExonicFunc, perl = TRUE), collapse = "|"),
      sf = paste0(SOMATICISM_FILTER, collapse = "|")
    ) %>%
    summarise(
      CHROM = first(CHROM),
      POS = first(POS),
      REF = first(REF),
      ALT = first(ALT),
      FILTER = first(FILTER),
      TRANSCRIPT = first(tr),
      AACHANGE = first(aa),
      GENEDETAIL = first(gd),
      FUNC = first(fc),
      EXONICFUNC = first(ef),
      FILTER.SOMATICISM = first(sf),
      vcf_idx = first(vcf_idx)
    )
  
  tmp_vcf_annot_vars_ann_som <- tmp_vcf_annot_vars_ann_som[c(
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "vcf_idx",
    "FILTER",
    "TRANSCRIPT",
    "AACHANGE",
    "GENEDETAIL",
    "FUNC",
    "EXONICFUNC",
    "FILTER.SOMATICISM"
  )]
  
  tmp_vcf_annot_vars_ann_som <- tmp_vcf_annot_vars_ann_som[
    order(
      tmp_vcf_annot_vars_ann_som$CHROM,
      tmp_vcf_annot_vars_ann_som$POS,
      tmp_vcf_annot_vars_ann_som$REF,
      tmp_vcf_annot_vars_ann_som$ALT,
      tmp_vcf_annot_vars_ann_som$vcf_idx
    ),
  ]
  
  tmp_vcf_annot_vars_ann_som <- tmp_vcf_annot_vars_ann_som %>%
    group_by(CHROM, POS, REF, ALT, FILTER) %>%
    mutate(
      sf = paste0(FILTER.SOMATICISM, collapse = ","),
      tr = paste0(TRANSCRIPT, collapse = ","),
      aa = paste0(AACHANGE, collapse = ","),
      gd = paste0(GENEDETAIL, collapse = ","),
      fc = paste0(FUNC, collapse = ","),
      ef = paste0(EXONICFUNC, collapse = ",")
    ) %>%
    summarise(
      FILTER.SOMATICISM = first(sf),
      INFO.ANN_SOM = paste(
        "ANNOTATION_METHOD=Annovar",
        paste("ANNOTATION_SOURCE", refseq_ensembl_suffix, sep = "="),
        paste("TRANSCRIPT", first(tr), sep = "="),
        paste("AACHANGE", first(aa), sep = "="),
        paste("GENEDETAIL", first(gd), sep = "="),
        paste("FUNC", first(fc), sep = "="),
        paste("EXONICFUNC", first(ef), sep = "="),
        sep = ";"
      )
    ) %>%
    mutate(
      INFO.SOMATICISM_FILTER = paste("SOMATICISM_FILTER", FILTER.SOMATICISM, sep = "=")
    )
  
  
  # ----- CHIP filter table ----- #
  tmp_vcf_annot_vars_chip_filtered_put <- vars_chip_filtered_put_df[c(
    "Chr",
    "Start",
    "End",
    "Ref",
    "Alt",
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "FILTER",
    "gene",
    transcript,
    aachange,
    genedetail,
    func,
    exonic_func,
    "mutation_class",
    "mutation_definition",
    "lineage",
    "putative",
    "chr",
    "strand",
    "c_term_genomic_start",
    "c_term_genomic_end",
    "publication_source_concat",
    "CHIP_FILTER",
    "PUTATIVE_FILTER",
    "COMBINED_FILTER"
  )]
  
  colnames(tmp_vcf_annot_vars_chip_filtered_put) <- c(
    "Chr",
    "Start",
    "End",
    "Ref",
    "Alt",
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "FILTER",
    "CHIP_gene",
    "CHIP_transcript",
    "CHIP_AAChange",
    "CHIP_GeneDetail",
    "CHIP_Func",
    "CHIP_ExonicFunc",
    "CHIP_mutation_class",
    "CHIP_mutation_definition",
    "CHIP_lineage",
    "CHIP_putative",
    "CHIP_gene_chr",
    "CHIP_gene_strand",
    "CHIP_gene_c_term_genomic_start",
    "CHIP_gene_c_term_genomic_end",
    "CHIP_publication_source_concat",
    "CHIP_FILTER",
    "PUTATIVE_FILTER",
    "COMBINED_FILTER"
  )
  
  # Only keep mutations that don't fail CHIP filter
  tmp_vcf_annot_vars_chip_filtered_put <- tmp_vcf_annot_vars_chip_filtered_put[tmp_vcf_annot_vars_chip_filtered_put$CHIP_FILTER != "FAIL",]
  
  tmp_vcf_annot_vars_chip_filtered_put$vcf_idx <- vars_annovar_vcf_ids_df[
    paste(
      tmp_vcf_annot_vars_chip_filtered_put$Chr,
      tmp_vcf_annot_vars_chip_filtered_put$Start,
      tmp_vcf_annot_vars_chip_filtered_put$End,
      tmp_vcf_annot_vars_chip_filtered_put$Ref,
      tmp_vcf_annot_vars_chip_filtered_put$Alt,
      sep = ":"
    ),
    "idx"
  ]
  
  tmp_vcf_annot_vars_chip_filtered_put <- tmp_vcf_annot_vars_chip_filtered_put %>%
    group_by(Chr, Start, End, Ref, Alt) %>%
    mutate(
      cg = paste0(CHIP_gene, collapse = "|"),
      ct = paste0(CHIP_transcript, collapse = "|"),
      ca = paste0(CHIP_AAChange, collapse = "|"),
      cgd = paste0(CHIP_GeneDetail, collapse = "|"),
      cfc = paste0(gsub("\\s+", "_", CHIP_Func, perl = TRUE), collapse = "|"),
      cef = paste0(gsub("\\s+", "_", CHIP_ExonicFunc, perl = TRUE), collapse = "|"),
      cc = paste0(CHIP_mutation_class, collapse = "|"),
      cd = paste0(CHIP_mutation_definition, collapse = "|"),
      cl = paste0(CHIP_lineage, collapse = "|"),
      cp = paste0(CHIP_putative, collapse = "|"),
      cgchr = paste0(CHIP_gene_chr, collapse = "|"),
      cgstr = paste0(CHIP_gene_strand, collapse = "|"),
      cgcs = paste0(CHIP_gene_c_term_genomic_start, collapse = "|"),
      cgce = paste0(CHIP_gene_c_term_genomic_end, collapse = "|"),
      cpub = paste0(gsub(";", "\\\\x3b", CHIP_publication_source_concat), collapse = "|"),
      cf = paste0(CHIP_FILTER, collapse = "|"),
      pf = paste0(PUTATIVE_FILTER, collapse = "|"),
      cmbf = paste0(COMBINED_FILTER, collapse = "|")
    ) %>%
    summarise(
      CHROM = first(CHROM),
      POS = first(POS),
      REF = first(REF),
      ALT = first(ALT),
      FILTER = first(FILTER),
      CHIP_gene = first(cg),
      CHIP_transcript = first(ct),
      CHIP_AAChange = first(ca),
      CHIP_GeneDetail = first(cgd),
      CHIP_Func = first(cfc),
      CHIP_ExonicFunc = first(cef),
      CHIP_mutation_class = first(cc),
      CHIP_mutation_definition = first(cd),
      CHIP_lineage = first(cl),
      CHIP_putative = first(cp),
      CHIP_gene_chr = first(cgchr),
      CHIP_gene_strand = first(cgstr),
      CHIP_gene_c_term_genomic_start = first(cgcs),
      CHIP_gene_c_term_genomic_end = first(cgce),
      CHIP_publication_source_concat = first(cpub),
      FILTER.CHIP = first(cf),
      FILTER.PUTATIVE = first(pf),
      FILTER.COMBINED = first(cmbf),
      vcf_idx = first(vcf_idx)
    )
  
  tmp_vcf_annot_vars_chip_filtered_put <- tmp_vcf_annot_vars_chip_filtered_put[c(
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "vcf_idx",
    "FILTER",
    "CHIP_gene",
    "CHIP_transcript",
    "CHIP_AAChange",
    "CHIP_GeneDetail",
    "CHIP_Func",
    "CHIP_ExonicFunc",
    "CHIP_mutation_class",
    "CHIP_mutation_definition",
    "CHIP_lineage",
    "CHIP_putative",
    "CHIP_gene_chr",
    "CHIP_gene_strand",
    "CHIP_gene_c_term_genomic_start",
    "CHIP_gene_c_term_genomic_end",
    "CHIP_publication_source_concat",
    "FILTER.CHIP",
    "FILTER.PUTATIVE",
    "FILTER.COMBINED"
  )]
  
  tmp_vcf_annot_vars_chip_filtered_put <- tmp_vcf_annot_vars_chip_filtered_put[
    order(
      tmp_vcf_annot_vars_chip_filtered_put$CHROM,
      tmp_vcf_annot_vars_chip_filtered_put$POS,
      tmp_vcf_annot_vars_chip_filtered_put$REF,
      tmp_vcf_annot_vars_chip_filtered_put$ALT,
      tmp_vcf_annot_vars_chip_filtered_put$vcf_idx
    ),
  ]
  
  tmp_vcf_annot_vars_chip_filtered_put <- tmp_vcf_annot_vars_chip_filtered_put %>%
    group_by(CHROM, POS, REF, ALT, FILTER) %>%
    mutate(
      cf = paste0(FILTER.CHIP, collapse = ","),
      pf = paste0(FILTER.PUTATIVE, collapse = ","),
      cmbf = paste0(FILTER.COMBINED, collapse = ","),
      cg = paste0(CHIP_gene, collapse = ","),
      ct = paste0(CHIP_transcript, collapse = ","),
      ca = paste0(CHIP_AAChange, collapse = ","),
      cgd = paste0(CHIP_GeneDetail, collapse = ","),
      cfc = paste0(CHIP_Func, collapse = ","),
      cef = paste0(CHIP_ExonicFunc, collapse = ","),
      cc = paste0(CHIP_mutation_class, collapse = ","),
      cd = paste0(CHIP_mutation_definition, collapse = ","),
      cl = paste0(CHIP_lineage, collapse = ","),
      cp = paste0(CHIP_putative, collapse = ","),
      cgchr = paste0(CHIP_gene_chr, collapse = ","),
      cgstr = paste0(CHIP_gene_strand, collapse = ","),
      cgcs = paste0(CHIP_gene_c_term_genomic_start, collapse = ","),
      cgce = paste0(CHIP_gene_c_term_genomic_end, collapse = ","),
      cpub = paste0(CHIP_publication_source_concat, collapse = ",")
    ) %>%
    summarise(
      FILTER.CHIP = first(cf),
      FILTER.PUTATIVE = first(pf),
      FILTER.COMBINED = first(cmbf),
      INFO.CHIP = paste(
        paste("CHIP_gene", first(cg), sep = "="),
        paste("CHIP_transcript", first(ct), sep = "="),
        paste("CHIP_AAChange", first(ca), sep = "="),
        paste("CHIP_GeneDetail", first(cgd), sep = "="),
        paste("CHIP_Func", first(cfc), sep = "="),
        paste("CHIP_ExonicFunc", first(cef), sep = "="),
        paste("CHIP_mutation_class", first(cc), sep = "="),
        paste("CHIP_mutation_definition", first(cd), sep = "="),
        paste("CHIP_lineage", first(cl), sep = "="),
        paste("CHIP_putative", first(cp), sep = "="),
        paste("CHIP_gene_chr", first(cgchr), sep = "="),
        paste("CHIP_gene_strand", first(cgstr), sep = "="),
        paste("CHIP_gene_c_term_genomic_start", first(cgcs), sep = "="),
        paste("CHIP_gene_c_term_genomic_end", first(cgce), sep = "="),
        paste("CHIP_publication_source_concat", first(cpub), sep = "="),
        sep = ";"
      )
    ) %>%
    mutate(
      INFO.CHIP_FILTER = paste("CHIP_FILTER", FILTER.CHIP, sep = "="),
      INFO.PUTATIVE_FILTER = paste("PUTATIVE_FILTER", FILTER.PUTATIVE, sep = "="),
      INFO.COMBINED_FILTER = paste("COMBINED_FILTER", FILTER.COMBINED, sep = "=")
    )
  
  vcf_annot_merge <- merge(vcf_annot, tmp_vcf_annot_vars_hf_hp, by = c("CHROM", "POS", "REF", "ALT", "FILTER"), all.x = TRUE)
  vcf_annot_merge <- merge(vcf_annot_merge, tmp_vcf_annot_vars_ann_som, by = c("CHROM", "POS", "REF", "ALT", "FILTER"), all.x = TRUE)
  vcf_annot_merge <- merge(vcf_annot_merge, tmp_vcf_annot_vars_chip_filtered_put, by = c("CHROM", "POS", "REF", "ALT", "FILTER"), all.x = TRUE)
  vcf_annot_merge <- as_tibble(vcf_annot_merge)
  
  vcf_annot_final <- vcf_annot_merge %>%
    mutate(
      INFO = gsub(";NA(?=;)", "", paste(
        INFO.start,
        INFO.HFHP,
        INFO.ANN_SOM,
        INFO.CHIP,
        INFO.HARD_FILTER,
        INFO.HOMOPOLYMER_FILTER,
        INFO.SOMATICISM_FILTER,
        INFO.CHIP_FILTER,
        INFO.PUTATIVE_FILTER,
        INFO.COMBINED_FILTER,
        INFO.end,
        sep = ";"
      ), perl = TRUE),
      ID = ".",
      QUAL = ".",
      FILTER.CHIP.PASS = if_else(
        grepl("PASS", FILTER.COMBINED, fixed = TRUE),
        "PASS",
        "FAIL",
        missing = "FAIL"
      ),
      FILTER.REDUCED.HARD = if_else(grepl("PASS", FILTER.HARD, fixed = TRUE), "", "CHIP_FAIL_HARD_FILTER", missing = ""),
      FILTER.REDUCED.HOMOPOLYMER = if_else(grepl("PASS", FILTER.HOMOPOLYMER, fixed = TRUE), "", "CHIP_FAIL_HOMOPOLYMER_FILTER", missing = ""),
      FILTER.REDUCED.SOMATICISM = if_else(grepl("PASS", FILTER.SOMATICISM, fixed = TRUE), "", "CHIP_FAIL_SOMATICISM_FILTER", missing = ""),
      FILTER.REDUCED.CHIP = if_else(grepl("PASS", FILTER.CHIP, fixed = TRUE), "", "CHIP_FAIL_CHIP_FILTER", missing = ""),
      FILTER.REDUCED.PUTATIVE = if_else(grepl("PASS", FILTER.PUTATIVE, fixed = TRUE), "", "CHIP_FAIL_PUTATIVE_FILTER", missing = "")
    ) %>%
    mutate(
      FILTER.CHIP.NOPASS = if_else(
        grepl("MANUAL_REVIEW", FILTER.COMBINED, fixed = TRUE) & FILTER.CHIP.PASS == "FAIL",
        "CHIP_MANUAL_REVIEW",
        gsub(
          "^;+", "", gsub(
            ";;+", ";", paste(
              FILTER.REDUCED.HARD,
              FILTER.REDUCED.HOMOPOLYMER,
              FILTER.REDUCED.SOMATICISM,
              FILTER.REDUCED.CHIP,
              FILTER.REDUCED.PUTATIVE,
              sep = ";"
            )
          )
        ),
        missing = "CHIP_FAIL"
      )
    ) %>%
    mutate(
      FILTER = if_else(
        FILTER.CHIP.PASS == "PASS",
        "PASS",
        gsub("^;+", "", gsub(";;+", ";", paste(gsub("^(\\.|PASS)$", "", FILTER), FILTER.CHIP.NOPASS, sep = ";"))),
        missing = gsub("^;+", "", gsub(";;+", ";", paste(gsub("^(\\.|PASS)$", "", FILTER), "CHIP_FAIL", sep = ";")))
      )
    )
  return(vcf_annot_final %>% select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO))
}
