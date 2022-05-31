# ============================================= #
# Define regex strings for mutation definitions #
# ============================================= #

chip_def_regex_list <- list(
  all = "^all$",
  c_terminal = "^c_term$",
  exon_any = "^exon\\d+$",
  exon_aa_insertion_specific = "^exon\\d+\\[ins[A-Z]\\]$",
  range = "^\\d+\\-\\d+$",
  substitution_specific = "^[A-Z]\\d+[A-Z]$",
  substitution_any = "^[A-Z]\\d+$",
  deletion_single = "^[A-Z]\\d+del$",
  deletion_range = "^[A-Z]\\d+_[A-Z]\\d+del$",
  insertion_specific = "^[A-Z]\\d+_[A-Z]\\d+ins[A-Z]+$",
  insertion_any = "^[A-Z]\\d+_[A-Z]\\d+ins$",
  deletion_insertion = "^[A-Z]\\d+_[A-Z]\\d+delins[A-Z]+$",
  coding_range = "^c\\.\\d+-\\d+$",
  coding_insertion_any = "^c\\.\\d+_\\d+ins$",
  frameshift_any = "^[A-Z]\\d+fs$",
  stop_gain = "^[A-Z]\\d+\\*$"
)

exonic_func_regex_list <- list(
  frameshift = "(^|;)frameshift ",
  stop_gain = "(^|;)stopgain(;|$)",
  nonsynonymous = "(^|;)(nonsynonymous|nonframeshift) "
)

gene_detail_regex <- list(
  sub = "^c\\.([\\-\\*])?(\\d+)([\\+\\-](\\d+))?([ACGT\\-]+)\\>([ACGT\\-]+)$",
  range = "^c\\.([\\-\\*])?(\\d+)([\\+\\-](\\d+))?_([\\-\\*])?(\\d+)([\\+\\-](\\d+))?(del|ins|delins)(\\-|[ACGT]+)$"
)

gene_detail_regex_groups <- list(
  sub = list(
    utr_intronic = 1,
    start = 2,
    end = 2,
    intron_offset_start = 4,
    intron_offset_end = 4,
    ref = 5,
    alt = 6,
    delins = NA
  ),
  range = list(
    utr_intronic = 1,
    start = 2,
    end = 6,
    intron_offset_start = 4,
    intron_offset_end = 8,
    ref = NA,
    alt = 10,
    delins = 9
  )
)

gene_detail_info_template <- list(
  utr_intronic = NA,
  start = NA,
  end = NA,
  intron_offset_start = NA,
  intron_offset_end = NA,
  ref = NA,
  alt = NA,
  delins = NA
)

aa_change_regex <- list(
  coding = list(
    sub = "^c\\.([ACGT])(\\d+)([ACGT])$",
    dupdelins = "^c\\.(\\d+)(dup|del|delins)([ACGT]+)$",
    range = "^c\\.(\\d+)_(\\d+)((del)|(dup|ins|delins)([ACGT]+))$"
  ),
  protein = list(
    sub = "^p\\.([A-Z\\*])(\\d+)([A-Z\\*])$",
    delins = "^p\\.([A-Z\\*])(\\d+)((del)|(delins)([A-Z\\*]+))$",
    delins_range = "^p\\.([A-Z\\*])(\\d+)_([A-Z\\*])?(\\d+)((del)|(ins|delins)([A-Z\\*]+))$",
    frameshift = "^p\\.([A-Z\\*])(\\d+)([A-Z\\*])fs\\*(\\d+|\\?)$",
    startloss = "^p\\.(M)(1)(\\?)$"
  )
)

aa_change_regex_groups <- list(
  coding = list(
    sub = list(
      start = 2,
      end = 2,
      dup = NA,
      del = NA,
      ins = NA,
      delins = NA,
      ref = 1,
      alt = 3
    ),
    dupdelins = list(
      start = 1,
      end = 1,
      dup = 2,
      del = 2,
      ins = NA,
      delins = 2,
      ref = 3,
      alt = 3  # this group is ref when del and alt when dup or delins
    ),
    range = list(
      start = 1,
      end = 2,
      dup = 5,
      del = 4,
      ins = 5,
      delins = 5,
      ref = NA,
      alt = 6  # no ref or alt present for del; alt present for dup, ins, and delins
    )
  ),
  protein = list(
    sub = list(
      start = 2,
      end = 2,
      ref_start = 1,
      ref_end = 1,
      alt = 3,
      del = NA,
      ins = NA,
      delins = NA,
      term = NA
    ),
    delins = list(
      start = 2,
      end = 2,
      ref_start = 1,
      ref_end = 1,
      alt = 6,
      del = 4,
      ins = NA,
      delins = 5,
      term = NA
    ),
    delins_range = list(
      start = 2,
      end = 4,
      ref_start = 1,
      ref_end = 3,
      alt = 8,
      del = 6,
      ins = 7,
      delins = 7,
      term = NA
    ),
    frameshift = list(
      start = 2,
      end = 2,
      ref_start = 1,
      ref_end = 1,
      alt = 3,
      del = NA,
      ins = NA,
      delins = NA,
      term = 4,
    ),
    startloss = list(
      start = 2,
      end = 2,
      ref_start = 1,
      ref_end = 1,
      alt = 3,
      del = NA,
      ins = NA,
      delins = NA,
      term = NA
    )
  )
)

aa_change_info_template <- list(
  coding <- list(
    start = NA,
    end = NA,
    dup = NA,
    del = NA,
    ins = NA,
    delins = NA,
    ref = NA,
    alt = NA
  ),
  protein <- list(
    start = NA,
    end = NA,
    ref_start = NA,
    ref_end = NA,
    alt = NA,
    del = NA,
    ins = NA,
    delins = NA,
    term = NA
  )
)


# ======================== #
# Main function definition #
# ======================== #

match_mut_def <- function(df, refseq_ensembl_suffix = "ensGene") {
  tmp_df <- df
  tmp_df$CHIP_FILTER <- NA
  
  # ===== Placeholder code =====
  ret <- tmp_df
  ret$CHIP_FILTER <- "PASS"  # Possible values: PASS, MANUAL_REVIEW, and FAIL
  # ===== End placeholder code ===== #
  
  # Update combined filter column
  ret$COMBINED_FILTER <- apply(ret[c("COMBINED_FILTER", "CHIP_FILTER")], 1, function(x) {
    if (x[[2]] == "PASS") {
      return(x[[1]])
    } else if (x[[2]] == "FAIL") {
      return("FAIL")
    } else if (x[[1]] != "FAIL") {
      # CHIP_FILTER can only be 'MANUAL_REVIEW' at this point
      # so if COMBINED_FILTER is not 'FAIL', return 'MANUAL_REVIEW', else return 'FAIL'
      return("MANUAL_REVIEW")
    } else {
      return("FAIL")
    }
  })
  return(ret)
}


# ==================== #
# Function definitions #
# ==================== #

parse_aa_change_mutation <- function(mutation, coding = FALSE) {
  if (coding) {
    # ===== Placeholder code =====
    return(aa_change_info_template$coding)
    # ===== End placeholder code ===== #
  } else {
    # ===== Placeholder code =====
    return(aa_change_info_template$protein)
    # ===== End placeholder code ===== #
  }
}

parse_aa_change <- function(aa_change) {
  info <- strsplit(aa_change, ":", fixed = TRUE)[[1]]
  info_list <- list(transcript = NA, exon = NA, mutation_coding = NA, mutation_protein = NA)
  info_list$transcript = gsub("\\.\\d+$", "", grep("^ENST\\d+", info, perl = TRUE, value = TRUE), perl = TRUE)
  info_list$exon <- gsub("^exon(\\d+)$", "\\1", grep("^exon", info, perl = TRUE, value = TRUE), perl = TRUE)
  info_list$mutation_coding <- grep("^c\\.", info, perl = TRUE, value = TRUE)
  info_list$mutation_protein <- grep("^p\\.", info, perl = TRUE, value = TRUE)
  info_list$mutation_coding_info <- parse_aa_change_mutation(info_list$mutation_coding, coding = TRUE)
  info_list$mutation_protein_info <- parse_aa_change_mutation(info_list$mutation_protein, coding = FALSE)
  return(info_list)
}

parse_gene_detail_mutation <- function(mutation) {
  # ===== Placeholder code =====
  return(gene_detail_info_template)
  # ===== End placeholder code ===== #
}

parse_gene_detail <- function(gene_detail) {
  info <- strsplit(gene_detail, ":", fixed = TRUE)[[1]]
  info_list <- list(transcript = NA, exon = NA, mutation = NA)
  info_list$transcript <- gsub("\\.\\d+$", "", grep("^ENST\\d+", info, perl = TRUE, value = TRUE), perl = TRUE)
  info_list$exon <- gsub("^exon(\\d+)$", "\\1", grep("^exon", info, perl = TRUE, value = TRUE), perl = TRUE)
  info_list$mutation <- grep("^c\\.", info, perl = TRUE, value = TRUE)
  return(info_list)
}

parse_chip_def <- function(chip_def) {
  
}



match_splicing <- function(df, refseq_ensembl_suffix = "ensGene") {
  tmp_df <- df
  if(!("CHIP_FILTER" %in% colnames(tmp_df))) {
    tmp_df$CHIP_FILTER <- NA
  }
  
  func <- paste("Func", refseq_ensembl_suffix, sep = ".")
  aachange <- paste("AAChange", refseq_ensembl_suffix, sep = ".")
  genedetail <- paste("GeneDetail", refseq_ensembl_suffix, sep = ".")
  
  splicing_idx <- tmp_df[[func]] == "splicing" & tmp_df$mutation_class == "splice_site"
  tmp_df_splicing <- tmp_df[splicing_idx,]
  tmp_df_nonsplicing <- tmp_df[!splicing_idx,]
  
  # ===== Do some stuff ===== #
  
  # ========================= #
  
  ret <- rbind(tmp_df_nonsplicing, tmp_df_splicing)
  return(ret)
}

match_exonic <- function(df, refseq_ensembl_suffix = "ensGene") {
  tmp_df <- df
  if(!("CHIP_FILTER" %in% colnames(tmp_df))) {
    tmp_df$CHIP_FILTER <- NA
  }
  
  func <- paste("Func", refseq_ensembl_suffix, sep = ".")
  aachange <- paste("AAChange", refseq_ensembl_suffix, sep = ".")

  exonic_idx <- tmp_df[[func]] == "exonic" & tmp_df$mutation_class != "splice_site"
  tmp_df_exonic <- tmp_df[exonic_idx,]
  tmp_df_nonexonic <- tmp_df[!exonic_idx,]
  
  # ===== Do some stuff ===== #
  
  # ========================= #
  
  ret <- rbind(tmp_df_nonexonic, tmp_df_exonic)
  return(ret)
}

# match_stop_gain <- function(df, refseq_ensembl_suffix = "ensGene") {
#   
# }
# 
# match_frameshift <- function(df, refseq_ensembl_suffix = "ensGene") {
#   
# }
# 
# match_nonsynonymous <- function(df, refseq_ensembl_suffix = "ensGene") {
#   
# }
