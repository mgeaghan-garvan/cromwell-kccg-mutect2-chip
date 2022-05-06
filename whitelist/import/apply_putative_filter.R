apply_putative_filter <- function(df, whitelist) {
  df$PUTATIVE <-  (whitelist$overall_wl$Putative_Whitelist | whitelist$overall_wl$Putative_Manual_Review)
  df$KNOWN <- (whitelist$overall_wl$Whitelist | whitelist$overall_wl$Manual_Review)
  df$PUTATIVE_FILTER <- apply(df[, c("PUTATIVE", "AD", "F1R2", "F2R1", "VAF", "KNOWN")], 1, function(x) {
    p <- as.logical(x[[1]])
    ad <- as.integer(strsplit(x[[2]], ",")[[1]])
    f1r2 <- as.integer(strsplit(x[[3]], ",")[[1]])
    f2r1 <- as.integer(strsplit(x[[4]], ",")[[1]])
    vaf <- as.numeric(strsplit(x[[5]], ",")[[1]])
    kwn <- as.logical(x[[6]])
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
    if (kwn) {
      PUTATIVE_FILTER = "MANUAL_REVIEW"
    } else {
      PUTATIVE_FILTER = "FAIL"
    }
    if(filter_pass) {
      PUTATIVE_FILTER = "PASS"
    } else if(filter_manual_review) {
      PUTATIVE_FILTER = "MANUAL_REVIEW"
    }
    return(PUTATIVE_FILTER)
  })
  df$COMBINED_FILTER <- apply(df[, c("COMBINED_FILTER", "PUTATIVE_FILTER")], 1, function(x) {
    if (x[[1]] == "PASS" && x[[2]] == "PASS") {
      return("PASS")
    } else if (x[[1]] == "FAIL" || x[[2]] == "FAIL") {
      return("FAIL")
    } else {
      return("MANUAL_REVIEW")
    }
  })
  return(df)
}