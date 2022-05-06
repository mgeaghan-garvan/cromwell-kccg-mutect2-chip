# ========================================================== #
# Define functions for parsing condition and AA change codes #
# ========================================================== #

parse_cond <- function(x) {
  cond_type <- NA
  start_pos <- NA
  end_pos <- NA
  AA_ref <- NA
  AA_alt <- NA
  exon_num <- NA
  if(grepl("^AA\\d+(\\-\\d+)?$", x, perl = TRUE)) {
    # Range of amino acid positions
    y <- str_match_all(x, "^AA(\\d+)(\\-(\\d+))?$")[[1]]
    start_pos <- as.integer(y[2])
    end_pos <- as.integer(y[4])
    if(is.na(end_pos)) {
      end_pos <- start_pos
    }
    cond_type <- "aa_range"
  } else if(grepl("^[A-Z]\\d+[A-Z]$", x, perl = TRUE)) {
    # AA substitution
    y <- str_match_all(x, "^([A-Z])(\\d+)([A-Z])$")[[1]]
    start_pos <- as.integer(y[3])
    end_pos <- start_pos
    AA_ref <- y[2]
    AA_alt <- y[4]
    cond_type <- "aa_change"
  } else if(grepl("^exon", x, perl = TRUE)) {
    # Exon number
    y <- str_match_all(x, "^exon(\\d+)$")[[1]]
    exon_num <- as.integer(y[2])
    cond_type <- "exon"
  }
  return(list(
    cond_type = cond_type,
    start_pos = start_pos,
    end_pos = end_pos,
    AA_ref = AA_ref,
    AA_alt = AA_alt,
    exon_num = exon_num
  ))
}

parse_aa_change <- function(x) {
  mut_type <- NA
  start_pos <- NA
  end_pos <- NA
  start_AA_ref <- NA
  start_AA_alt <- NA
  end_AA_ref <- NA
  if(grepl("^[A-Z]\\d+[A-Z\\*]$", x, perl = TRUE)) {
    # AA substitution
    y <- str_match_all(x, "^([A-Z])(\\d+)([A-Z\\*])$")[[1]]
    start_pos <- as.integer(y[3])
    end_pos <- start_pos
    start_AA_ref <- y[2]
    start_AA_alt <- y[4]
    mut_type <- "aa_change"
  } else if(grepl("^[A-Z]\\d+(_[A-Z]\\d+)?(ins|del)",x, perl = TRUE)) {
    # Non-frameshift INDEL
    y <- str_match_all(x, "^([A-Z])(\\d+)(_([A-Z])(\\d+))?(ins|del)")[[1]]
    start_pos <- as.integer(y[3])
    end_pos <- as.integer(y[6])
    if(is.na(end_pos)) {
      end_pos <- start_pos
    }
    start_AA_ref <- y[2]
    end_AA_ref <- y[5]
    mut_type <- y[7]
  } else if(grepl("^[A-Z]\\d+[A-Z]fs.*$", x, perl = TRUE)) {
    # Frameshift mutation
    y <- str_match_all(x, "^([A-Z])(\\d+)([A-Z])(fs).*$")[[1]]
    start_pos <- as.integer(y[3])
    end_pos <- start_pos
    start_AA_ref <- y[2]
    start_AA_alt <- y[4]
    mut_type <- y[5]
  }
  return(list(
    mut_type = mut_type,
    start_pos = start_pos,
    end_pos = end_pos,
    start_AA_ref = start_AA_ref,
    start_AA_alt = start_AA_alt,
    end_AA_ref = end_AA_ref
  ))
}