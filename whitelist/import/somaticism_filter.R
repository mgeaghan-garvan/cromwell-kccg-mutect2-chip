apply_somaticism_filter <- function(df, somaticism_filter_file, transcript_col) {
  tmp_df <- df
  
  # Read in somaticism file
  somaticism_list <- scan(somaticism_file, character())
  
  som_select <- tmp_df[[transcript_col]] %in% somaticism_list
  
  tmp_df_som <- tmp_df[som_select,]
  tmp_df_non_som <- tmp_df[!som_select,]
  
  if (dim(tmp_df_som)[1] == 0) {
    tmp_df$SOMATICISM_FILTER <- "PASS"
    return(tmp_df)
  }
  
  tmp_df_som$SOMATICISM_FILTER <- unlist(apply(tmp_df_som[c("AD", "DP")], 1, function(x) {
    ad <- as.integer(strsplit(x[[1]], ",")[[1]])
    ad <- ad[2:length(ad)]
    dp <- as.integer(x[[2]])
    
    bt <- unlist(lapply(ad, function(x) { binom.test(x, dp, 0.5)$p.value }))
    
    filter_pass <- all(bt < 0.001)
    filter_manual_review <- any(bt < 0.001)
    
    SOMATICISM_FILTER = "FAIL"
    
    if(filter_pass) {
      SOMATICISM_FILTER = "PASS"
    } else if(filter_manual_review) {
      SOMATICISM_FILTER = "MANUAL_REVIEW"
    }
    
    return(SOMATICISM_FILTER)
  }))
  
  if (dim(tmp_df_non_som)[1] == 0) {
    tmp_df_non_som$SOMATICISM_FILTER <- character(0)
  } else {
    tmp_df_non_som$SOMATICISM_FILTER <- "PASS"
  }
  
  tmp_df_2 <- rbind(tmp_df_som, tmp_df_non_som)
  
  # Update combined filter column
  tmp_df_2$COMBINED_FILTER <- apply(tmp_df_2[c("COMBINED_FILTER", "SOMATICISM_FILTER")], 1, function(x) {
    if (x[[2]] == "PASS") {
      return(x[[1]])
    } else if (x[[2]] == "FAIL") {
      return("FAIL")
    } else if (x[[1]] != "FAIL") {
      return("MANUAL_REVIEW")
    } else {
      return("FAIL")
    }
  })
  
  return(tmp_df_2)
}
