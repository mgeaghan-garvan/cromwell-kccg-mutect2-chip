library(stringr)

# ============================================= #
# Define regex strings for mutation definitions #
# ============================================= #

mutation_info_template = list(
  start = NA,
  end = NA,
  ref_start = NA,
  ref_end = NA,
  alt = NA,
  exon = NA,
  dupdelins = NA,
  utr_intronic = NA,
  intron_offset_start = NA,
  intron_offset_end = NA,
  fs_term = NA
)

chip_def_regex_list <- list(
  all = "^all$",
  c_terminal = "^c_term$",
  exon = "^exon(\\d+)$",
  exon_aa_insertion = "^exon(\\d+)\\[(ins)([A-Z])\\]$",
  range = "^(\\d+)\\-(\\d+)$",
  substitution = "^([A-Z])(\\d+)([A-Z])?$",
  deletion = "^([A-Z])(\\d+)(_([A-Z])(\\d+))?(del)$",
  insertion = "^([A-Z])(\\d+)_([A-Z])(\\d+)(ins)([A-Z]+)?$",
  deletion_insertion = "^([A-Z])(\\d+)(_([A-Z])(\\d+))?(delins)([A-Z]+)$",
  coding_range = "^c\\.(\\d+)-(\\d+)$",
  coding_insertion = "^c\\.(\\d+)_(\\d+)(ins)$",
  frameshift = "^([A-Z])(\\d+)fs$",
  stop_gain = "^([A-Z])(\\d+)\\*$"
)

chip_def_regex_groups <- list(
  all = list(),
  c_terminal = list(),
  exon = list(
    exon = 1
  ),
  exon_aa_insertion = list(
    alt = 3,
    exon = 1,
    dupdelins = 2
  ),
  range = list(
    start = 1,
    end = 2
  ),
  substitution = list(
    start = 2,
    end = 2,
    ref_start = 1,
    ref_end = 1,
    alt = 3
  ),
  deletion = list(
    start = 2,
    end = c(5, 2),
    ref_start = 1,
    ref_end = c(4, 1),
    dupdelins = 6
  ),
  insertion = list(
    start = 2,
    end = 4,
    ref_start = 1,
    ref_end = 3,
    alt = 6,
    dupdelins = 5
  ),
  deletion_insertion = list(
    start = 2,
    end = c(5, 2),
    ref_start = 1,
    ref_end = c(4, 1),
    alt = 7,
    dupdelins = 6
  ),
  coding_range = list(
    start = 1,
    end = 2
  ),
  coding_insertion = list(
    start = 1,
    end = 2,
    dupdelins = 3
  ),
  frameshift = list(
    start = 2,
    end = 2,
    ref_start = 1,
    ref_end = 1
  ),
  stop_gain = list(
    start = 2,
    end = 2,
    ref_start = 1,
    ref_end = 1
  )
)

exonic_func_regex_list <- list(
  frameshift = "(^|;)frameshift ",
  stop_gain = "(^|;)stopgain(;|$)",
  nonsynonymous = "(^|;)(nonsynonymous|nonframeshift) "
)

gene_detail_regex_list <- list(
  sub = "^c\\.([\\-\\*])?(\\d+)([\\+\\-](\\d+))?([ACGT\\-]+)\\>([ACGT\\-]+)$",
  range = "^c\\.([\\*])?(\\d+)?([\\+\\-](\\d+))?_([\\*])?(\\d+)?([\\+\\-](\\d+))?(del|ins|delins)(\\-|[ACGT]+)$"
)

gene_detail_regex_groups <- list(
  sub = list(
    utr_intronic = 1,
    start = 2,
    end = 2,
    intron_offset_start = 4,
    intron_offset_end = 4,
    ref_start = 5,
    ref_end = 5,
    alt = 6
  ),
  range = list(
    utr_intronic = 1,
    start = 2,
    end = 6,
    intron_offset_start = 4,
    intron_offset_end = 8,
    alt = 10,
    delins = 9
  )
)

aa_change_regex_list <- list(
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
      ref_start = 1,
      ref_end = 1,
      alt = 3
    ),
    dupdelins = list(
      start = 1,
      end = 1,
      dupdelins = 2,
      ref_start = 3,
      ref_end = 3,
      alt = 3  # this group is ref when del and alt when dup or delins
    ),
    range = list(
      start = 1,
      end = 2,
      dupdelins = c(4, 5),
      alt = 6  # no ref or alt present for del; alt present for dup, ins, and delins
    )
  ),
  protein = list(
    sub = list(
      start = 2,
      end = 2,
      ref_start = 1,
      ref_end = 1,
      alt = 3
    ),
    delins = list(
      start = 2,
      end = 2,
      ref_start = 1,
      ref_end = 1,
      alt = 6,
      dupdelins = c(4, 5)
    ),
    delins_range = list(
      start = 2,
      end = 4,
      ref_start = 1,
      ref_end = 3,
      alt = 8,
      dupdelins = c(6, 7)
    ),
    frameshift = list(
      start = 2,
      end = 2,
      ref_start = 1,
      ref_end = 1,
      alt = 3,
      fs_term = 4
    ),
    startloss = list(
      start = 2,
      end = 2,
      ref_start = 1,
      ref_end = 1,
      alt = 3
    )
  )
)


# =========================== #
# Helper function definitions #
# =========================== #

parse_mutation <- function(mutation, regex_list, regex_groups, exon = NA) {
  ret <- mutation_info_template
  for (mut_type in names(regex_list)) {
    # Iterate through each (mutally exclusive) mutation regex string
    rgx <- regex_list[[mut_type]]
    # Try and match the mutation string to the current regex string
    m <- str_match(mutation, rgx)
    if (all(is.na(m))) {
      # Move to next regex string if no match
      next
    }
    rgx_groups <- regex_groups[[mut_type]]
    for (field in names(rgx_groups)) {
      # Iterate through the regex capture groups defined for the current mutation type
      idx <- rgx_groups[[field]]
      if (is.na(idx) || !is.numeric(idx)) {
        next
      }
      for (i in idx) {
        # Get regex capture group value
        value <- m[1, (i + 1)]  # 1st capture group is in the second column, etc.
        if (!is.na(value)) {
          ret[[field]] <- value
          # Multiple group indices for a single field should be mutually exclusive, so we can break the loop
          break
        }
      }
    }
    # Break, since regex strings are mutually exclusive
    break
  }
  if (!is.na(exon)) {
    ret$exon <- exon
  }
  return(ret)
}

parse_aa_change <- function(aa_change) {
  info <- strsplit(aa_change, ":", fixed = TRUE)[[1]]
  info_list <- list(transcript = NA, mutation_class = NA, mutation_coding = NA, mutation_protein = NA, nonsynonymous_type = NA)
  info_list$transcript = gsub("\\.\\d+$", "", grep("^ENST\\d+", info, perl = TRUE, value = TRUE), perl = TRUE)
  info_list$mutation_coding <- grep("^c\\.", info, perl = TRUE, value = TRUE)
  info_list$mutation_protein <- grep("^p\\.", info, perl = TRUE, value = TRUE)
  exon <- gsub("^exon(\\d+)$", "\\1", grep("^exon", info, perl = TRUE, value = TRUE), perl = TRUE)
  info_list$mutation_coding_info <- parse_mutation(info_list$mutation_coding, aa_change_regex_list$coding, aa_change_regex_groups$coding, exon = exon)
  info_list$mutation_protein_info <- parse_mutation(info_list$mutation_protein, aa_change_regex_list$protein, aa_change_regex_groups$protein, exon = exon)
  if (grepl(aa_change_regex$protein$frameshift, info_list$mutation_protein, perl = TRUE)) {
    info_list$mutation_class <- "frameshift"
  } else if (!grepl("[X\\*]$", info_list$mutation_protein_info$ref_end, perl = TRUE) &&
             grepl("[X\\*]$", info_list$mutation_protein_info$alt, perl = TRUE)) {
    info_list$mutation_class <- "stop_gain"
  } else if (!all(is.na(info_list$mutation_protein_info)) &&
             (
               is.na(info_list$mutation_protein_info$ref_start) ||
               is.na(info_list$mutation_protein_info$alt) ||
               info_list$mutation_protein_info$ref_start != info_list$mutation_protein_info$alt
             )) {
    info_list$mutation_class <- "nonsynonymous"
    if (grepl(aa_change_regex$protein$sub, info_list$mutation_protein, perl = TRUE)) {
      info_list$nonsynonymous_type <- "substitution"
    } else if (!is.na(info_list$mutation_protein_info$dupdelins)) {
      info_list$nonsynonymous_type <- info_list$mutation_protein_info$dupdelins
    }
  }
  return(info_list)
}

parse_gene_detail <- function(gene_detail) {
  info <- strsplit(gene_detail, ":", fixed = TRUE)[[1]]
  info_list <- list(transcript = NA, mutation = NA, mutation_class = NA)
  info_list$transcript <- gsub("\\.\\d+$", "", grep("^ENST\\d+", info, perl = TRUE, value = TRUE), perl = TRUE)
  info_list$mutation <- grep("^c\\.", info, perl = TRUE, value = TRUE)
  exon <- gsub("^exon(\\d+)$", "\\1", grep("^exon", info, perl = TRUE, value = TRUE), perl = TRUE)
  info_list$mutation_coding_info <- parse_mutation(info_list$mutation, gene_detail_regex_list, gene_detail_regex_groups, exon = exon)
  if (!all(is.na(info_list$mutation_coding_info))) {
    info_list$mutation_class <- "splice_site"
  }
  return(info_list)
}

parse_chip_def <- function(chip_def, mut_class = NA) {
  info_list <- list()
  info_list$mutation_class <- mut_class
  info_list$mutation <- chip_def
  info_list$mutation_type <- NA
  info_list$mutation_info <- mutation_info_template
  for (mut_type in names(chip_def_regex_list)) {
    # Iterate through each (mutally exclusive) chip definition regex string
    rgx <- chip_def_regex_list[[mut_type]]
    # Try and match the chip definition string to the current regex string
    m <- str_match(info_list$mutation, rgx)
    if (all(is.na(m))) {
      # Move to next regex string if no match
      next
    }
    # Record mutation type
    info_list$mutation_type <- mut_type
    rgx_groups <- chip_def_regex_groups[[mut_type]]
    for (field in names(rgx_groups)) {
      # Iterate through the regex capture groups defined for the current mutation type
      idx <- rgx_groups[[field]]
      if (is.na(idx) || !is.numeric(idx)) {
        next
      }
      for (i in idx) {
        # Get regex capture group value
        value <- m[1, (i + 1)]  # 1st capture group is in the second column, etc.
        if (!is.na(value)) {
          info_list$mutation_info[[field]] <- value
          # Multiple group indices for a single field should be mutually exclusive, so we can break the loop
          break
        }
      }
    }
    # Break, since regex strings are mutually exclusive
    break
  }
  return(info_list)
}


# ======================== #
# Main function definition #
# ======================== #

match_mut_def <- function(df, refseq_ensembl_suffix = "ensGene") {
  aachange <- paste("AAChange", refseq_ensembl_suffix, sep = ".")
  genedetail <- paste("GeneDetail", refseq_ensembl_suffix, sep = ".")
  func <- paste("Func", refseq_ensembl_suffix, sep = ".")

  tmp_df <- df
  tmp_df$CHIP_FILTER <- unlist(apply(tmp_df[c(genedetail, aachange, func, "Chr", "Start", "End", "mutation_class", "mutation_definition", "chr", "c_term_genomic_start", "c_term_genomic_end")], 1, function(x) {
    gene_detail <- as.character(x[[1]])
    aa_change <- as.character(x[[2]])
    mut_func <- as.character(x[[3]])
    mut_chr <- as.character(x[[4]])
    mut_start <- as.character(x[[5]])
    mut_end <- as.character(x[[6]])
    chip_mut_class <- as.character(x[[7]])
    chip_mut_def <- as.character(x[[8]])
    chip_mut_chr <- as.character(x[[9]])
    chip_mut_c_term_start <- as.numeric(x[[10]])
    chip_mut_c_term_end <- as.numeric(x[[11]])
    
    gene_detail_info <- parse_gene_detail(gene_detail)
    aa_change_info <- parse_aa_change(aa_change)
    chip_info <- parse_chip_def(chip_mut_def, chip_mut_class)
    
    # Determine mutation class
    mut_class <- NA
    if (mut_func == "splicing") {
      mut_class <- "splice_site"
    } else {
      mut_class <- aa_change_info$mutation_class
    }
    
    # Handle unclassified mutations
    if (is.na(mut_class)) {
      return("FAIL")
    }
    
    # Handle nonsynonymous missense mutations
    if (chip_mut_class == "missense" &&
        mut_class == "nonsynonymous" &&
        !is.na(aa_change_info$nonsynonymous_type) &&
        aa_change_info$nonsynonymous_type == "substitution" &&
        !is.na(aa_change_info$mutation_protein_info$ref_start) &&
        !is.na(aa_change_info$mutation_protein_info$alt)) {
      mut_class <- "missense"
    }
    
    # Handle mismatched mutation class
    if (mut_class != chip_mut_class) {
      return("FAIL")
    }
    
    # Handle C-terminal mutations - do this FIRST
    mut_in_c_term <- FALSE
    if (
      mut_chr == chip_mut_chr && (
        (chip_mut_c_term_start <= mut_start && mut_start <= chip_mut_c_term_end) ||
        (chip_mut_c_term_start <= mut_end && mut_end <= chip_mut_c_term_end)
      )
    ) {
      mut_in_c_term <- TRUE
    }
    if (chip_mut_def == "c_term") {
      if (mut_in_c_term) {
        # C-terminal mutations pass if the CHIP mutation is defined as a C-terminal mutation
        return("PASS")
      } else {
        # Non-C-terminal mutations fail if the CHIP mutation is defined as a C-terminal mutation
        return("FAIL")
      }
    } else if (mut_in_c_term) {
      # C-terminal mutations fail if the CHIP mutation is NOT defined as a C-terminal mutation
      return("FAIL")
    }
    
    # Handle "all" CHIP mutation definitions
    if (chip_mut_def == "all") {
      return("PASS")
    }
    
    # Handle CHIP mutations defined by any mutation in a given exon
    if (chip_info$mutation_type == "exon") {
      if (mut_class == "splice_site" && !is.na(gene_detail_info$mutation_coding_info$exon && gene_detail_info$mutation_coding_info$exon == chip_info$mutation_info$exon)) {
        return("PASS")
      } else if (!is.na(aa_change_info$mutation_protein_info$exon) && aa_change_info$mutation_protein_info$exon == chip_info$mutation_info$exon) {
        return("PASS")
      } else {
        return("FAIL")
      }
    }
    
    # Handle CHIP mutations defined by specific amino acid insertions in a given exon
    if (chip_info$mutation_type == "exon_aa_insertion") {
      if (!is.na(aa_change_info$mutation_protein_info$exon) && aa_change_info$mutation_protein_info$exon == chip_info$mutation_info$exon) {
        # Mutation has valid exon recorded and matches the CHIP definition exon
        if (!is.na(aa_change_info$mutation_protein_info$dupdelins) && aa_change_info$mutation_protein_info$dupdelins == "ins" && !is.na(aa_change_info$mutation_protein_info$alt) && aa_change_info$mutation_protein_info$alt == chip_info$mutation_info$alt) {
          # Mutation is an insertion and the alternate allele matches the CHIP definition
          return("PASS")
        }
      }
      return("FAIL")
    }
    
    # Handle amino acid ranges
    if (chip_info$mutation_type == "range") {
      if (!is.na(aa_change_info$mutation_protein_info$start) &&
          !is.na(aa_change_info$mutation_protein_info$end) &&
          (
            (chip_info$mutation_info$start <= aa_change_info$mutation_protein_info$start && aa_change_info$mutation_protein_info$start <= chip_info$mutation_info$end) ||
            (chip_info$mutation_info$start <= aa_change_info$mutation_protein_info$end && aa_change_info$mutation_protein_info$end <= chip_info$mutation_info$end)
          )) {
        # Passes if mutation amino acid range overlaps with defined CHIP mutation amino acid range
        return("PASS")
      }
    }
    
    # Handle simple substitutions
    if (chip_info$mutation_type == "substitution" && !is.na(aa_change_info$nonsynonymous_type) && aa_change_info$nonsynonymous_type == "sub") {
      
    }
    
    # Handle deletions
    if (chip_info$mutation_type == "deletion" && !is.na(aa_change_info$nonsynonymous_type) && aa_change_info$nonsynonymous_type == "del") {
      
    }
    
    # Handle insertions
    if (chip_info$mutation_type == "insertion" && !is.na(aa_change_info$nonsynonymous_type) && aa_change_info$nonsynonymous_type == "ins") {
      
    }
    
    # Handle deletion-insertion events
    if (chip_info$mutation_type == "deletion_insertion" && !is.na(aa_change_info$nonsynonymous_type) && aa_change_info$nonsynonymous_type == "delins") {
      
    }
    
    # Handle coding sequence ranges
    if (chip_info$mutation_type == "coding_range") {
      if (mut_class == "splice_site" && !all(is.na(gene_detail_info$mutation_coding_info))) {
        
      } else if (!all(is.na(aa_change_info$mutation_coding_info))) {
        
      } else {
        return("FAIL")
      }
    }
    
    # Handle coding insertions
    if (chip_info$mutation_type == "coding_insertion") {
      if (mut_class == "splice_site" && !is.na(gene_detail_info$mutation_coding_info$dupdelins) && gene_detail_info$mutation_coding_info$dupdelins == "ins") {
        
      } else if (!is.na(aa_change_info$mutation_coding_info) && aa_change_info$mutation_coding_info == "ins") {
        
      } else {
        return("FAIL")
      }
    }
    
    # Handle frameshift events
    if (chip_info$mutation_type == "frameshift" && mut_class == "frameshift") {
      
    }
    
    # Handle stop gain events
    if (chip_info$mutation_type == "stop_gain" && mut_class == "stop_gain") {
      
    }
    
    # If there haven't been any matches at this point, default to FAIL
    return("FAIL")
  }))
  
  # ===== Placeholder code =====
  tmp_df$CHIP_FILTER <- "PASS"  # Possible values: PASS, MANUAL_REVIEW, and FAIL
  # ===== End placeholder code ===== #
  
  # Update combined filter column
  tmp_df$COMBINED_FILTER <- apply(tmp_df[c("COMBINED_FILTER", "CHIP_FILTER")], 1, function(x) {
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
  return(tmp_df)
}
