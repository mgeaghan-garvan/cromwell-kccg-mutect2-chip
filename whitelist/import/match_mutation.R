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
  frameshift = "^frameshift ",
  stop_gain = "^stopgain$",
  nonsynonymous = "^(nonsynonymous|nonframeshift) "
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

match_splicing <- function(df, refseq_ensembl_suffix = "ensGene") {
  tmp_df <- df
  if(!("CHIP_FILTER" %in% colnames(tmp_df))) {
    tmp_df$CHIP_FILTER <- NA
  }
  
  func <- paste("Func", refseq_ensembl_suffix, sep = ".")
  aachange <- paste("AAChange", refseq_ensembl_suffix, sep = ".")
  genedetail <- paste("GeneDetail", refseq_ensembl_suffix, sep = ".")
  
  splicing_idx <- grepl("splicing", tmp_df[[func]], fixed = TRUE) & tmp_df$mutation_class == "splice_site"
  tmp_df_splicing <- tmp_df[splicing_idx,]
  tmp_df_nonsplicing <- tmp_df[!splicing_idx,]
  
  # ===== Do some stuff ===== #
  
  # ========================= #
  
  ret <- rbind(tmp_df_nonsplicing, tmp_df_splicing)
}
