apply_putative_filter <- function(df) {
  tmp_df <- df
  if (dim(tmp_df)[1] == 0) {
    tmp_df$PUTATIVE_FILTER <- character(0)
    return(tmp_df)
  }
  tmp_df$PUTATIVE_FILTER <- apply(tmp_df[, c("putative", "AD", "F1R2", "F2R1", "VAF")], 1, function(x) {
    p <- as.logical(x[[1]])
    ad <- as.integer(strsplit(x[[2]], ",")[[1]])
    f1r2 <- as.integer(strsplit(x[[3]], ",")[[1]])
    f2r1 <- as.integer(strsplit(x[[4]], ",")[[1]])
    vaf <- as.numeric(strsplit(x[[5]], ",")[[1]])

    if (!p) {
      return("PASS")
    }
    
    # Allelic depth
    var_n_ad_gte_5 <- ad >= 5
    # VAF
    var_n_vaf_lte_20pc <- vaf <= 0.2
    # Forward and reverse strand support
    var_n_for_rev_gte_2 <- all(f1r2 >= 2) && all(f2r1 >= 2)
    
    filter_pass <- (
      all(var_n_ad_gte_5[2:length(var_n_ad_gte_5)]) &  # first allele = ref
        all(var_n_vaf_lte_20pc) &
        var_n_for_rev_gte_2
    )
    filter_manual_review <- (
      any(var_n_ad_gte_5[2:length(var_n_ad_gte_5)]) &  # first allele = ref; if any, but not all pass, go to manual review
        any(var_n_vaf_lte_20pc) &  # some but not all are passing
        var_n_for_rev_gte_2
    )
    
    PUTATIVE_FILTER = "FAIL"
    
    if(filter_pass) {
      PUTATIVE_FILTER = "PASS"
    } else if(filter_manual_review) {
      PUTATIVE_FILTER = "MANUAL_REVIEW"
    }
    
    return(PUTATIVE_FILTER)
  })
  
  # Update combined filter column
  tmp_df$COMBINED_FILTER <- apply(tmp_df[c("COMBINED_FILTER", "PUTATIVE_FILTER")], 1, function(x) {
    if (x[[2]] == "PASS") {
      return(x[[1]])
    } else if (x[[2]] == "FAIL") {
      return("FAIL")
    } else if (x[[1]] != "FAIL") {
      # PUTATIVE_FILTER can only be 'MANUAL_REVIEW' at this point
      # so if COMBINED_FILTER is not 'FAIL', return 'MANUAL_REVIEW', else return 'FAIL'
      return("MANUAL_REVIEW")
    } else {
      return("FAIL")
    }
  })
  
  return(tmp_df)
}