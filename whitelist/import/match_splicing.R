# =========================================================== #
# Define functions for matching variants with CHIP conditions #
# =========================================================== #

match_sp_conditions <- function(conditions, transcript_changes, transcript_exons, parse_cond_func) {
  whitelist <- FALSE
  manual_review <- FALSE

    for(cond in conditions) {
      if(cond == "any") {
        whitelist <- TRUE
        break
      } else if(cond == "c_term") {
        # No easy way to determine whether variant is in the C-terminal region - manual review
        manual_review <- TRUE
        next
      } else if(cond == "nan" || cond == "" || is.na(cond)) {
        # No nonsynonymous condition defined for this gene
        next
      } else if(transcript_changes == "nan" || transcript_changes == "" || is.na(transcript_changes) || transcript_exons == "nan" || transcript_exons == "" || is.na(transcript_exons)) {
        # Variant doesn't match the representative transcript(s) - cannot match to specific change or exon
        next
      }
      cond_details <- parse_cond_func(cond)
      if(all(is.na(cond_details))) {
        next
      } else if(cond_details$cond_type == "exon" && !is.na(cond_details$exon_num) && cond_details$exon_num %in% transcript_exons) { #####
        # Matching nonsynonymous variant in exon
        whitelist <- TRUE
        break
      } else if(cond_details$cond_type %in% c("aa_change", "aa_range")) {
        # AA substitutions and ranges aren't defined for splicing variants - manual review
        manual_review <- TRUE
        next
      } else if(is.na(cond_details$start_pos) || is.na(cond_details$end_pos)) {
        next
      }
    }
    if(whitelist) {
      manual_review <- FALSE
    }

  return(list(whitelist = whitelist, manual_review = manual_review))
}

match_sp <- function(x, parse_cond_func) {
  whitelist <- FALSE
  manual_review <- FALSE
  putative_whitelist <- FALSE
  putative_manual_review <- FALSE
  if(grepl("splicing", x[[1]], perl = TRUE)) {
    conditions <- strsplit(x[[2]], ";")[[1]]
    conditions_put <- strsplit(x[[3]], ";")[[1]]
    transcript_changes <- strsplit(x[[4]], ";")[[1]]
    transcript_exons <- as.integer(strsplit(x[[5]], ";")[[1]])
    match_conditions <- match_sp_conditions(conditions, transcript_changes, transcript_exons, parse_cond_func)
    match_put_conditions <- match_sp_conditions(conditions_put, transcript_changes, transcript_exons, parse_cond_func)
    whitelist <- match_conditions$whitelist
    manual_review <- match_conditions$manual_review
    putative_whitelist <- match_put_conditions$whitelist
    putative_manual_review <- match_put_conditions$manual_review
  }
  return(c(whitelist, manual_review, putative_whitelist, putative_manual_review))
}