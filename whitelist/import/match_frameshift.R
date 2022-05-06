# =========================================================== #
# Define functions for matching variants with CHIP conditions #
# =========================================================== #

match_fs_conditions <- function(conditions, aa_changes, aa_exons, parse_cond_func, parse_aa_change_func) {
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
    } else if(aa_changes == "nan" || aa_changes == "" || is.na(aa_changes) || aa_exons == "nan" || aa_exons == "" || is.na(aa_exons)) {
      # Variant doesn't match the representative transcript(s) - cannot match to specific AA or exon
      next
    }
    cond_details <- parse_cond_func(cond)
    if(all(is.na(cond_details))) {
      next
    } else if(cond_details$cond_type == "exon" && !is.na(cond_details$exon_num) && cond_details$exon_num %in% aa_exons) {
      # Matching nonsynonymous variant in exon
      whitelist <- TRUE
      break
    } else if(is.na(cond_details$start_pos) || is.na(cond_details$end_pos)) {
      next
    }
    for(a in aa_changes) {
      a_details <- parse_aa_change_func(a)
      if(all(is.na(a_details))) {
        next
      }
      if(a_details$mut_type == "fs" & (a_details$start_pos %in% cond_details$start_pos:cond_details$end_pos)) {
        # Current variant is a frameshift mutation
        if(cond_details$cond_type == "aa_range") {
          # Current variant in defined range
          whitelist <- TRUE
          break
        } else if(cond_details$cond_type == "aa_change" && !is.na(cond_details$AA_ref) && !is.na(cond_details$AA_alt)) {
          if(a_details$start_AA_ref == cond_details$AA_ref && a_details$start_AA_alt == cond_details$AA_alt) {
            # Current variant matches defined variant
            whitelist <- TRUE
            break
          }
        }
      }
    }
    if(whitelist) {
      break
    }
  }
  if(whitelist) {
    manual_review <- FALSE
  }
  return(list(whitelist = whitelist, manual_review = manual_review))
}

match_fs <- function(x, parse_cond_func, parse_aa_change_func) {
  whitelist <- FALSE
  manual_review <- FALSE
  putative_whitelist <- FALSE
  putative_manual_review <- FALSE
  if(grepl("^frameshift", x[[1]], perl = TRUE)) {
    conditions <- strsplit(x[[2]], ";")[[1]]
    conditions_put <- strsplit(x[[3]], ";")[[1]]
    aa_changes <- strsplit(x[[4]], ",")[[1]]
    aa_exons <- as.integer(strsplit(x[[5]], ",")[[1]])
    match_conditions <- match_fs_conditions(conditions, aa_changes, aa_exons, parse_cond_func, parse_aa_change_func)
    match_put_conditions <- match_fs_conditions(conditions_put, aa_changes, aa_exons, parse_cond_func, parse_aa_change_func)
    whitelist <- match_conditions$whitelist
    manual_review <- match_conditions$manual_review
    putative_whitelist <- match_put_conditions$whitelist
    putative_manual_review <- match_put_conditions$manual_review
  }
  return(c(whitelist, manual_review, putative_whitelist, putative_manual_review))
}