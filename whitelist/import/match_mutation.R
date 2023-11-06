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
  fs_term = NA,
  mutation_pattern = NA
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
  stop_gain_substitution = "^([A-Z])(\\d+)\\*$"
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
  stop_gain_substitution = list(
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
  sub = "^c\\.([\\*])?(\\d+)?([\\+\\-](\\d+))?([ACGT\\-]+)\\>([ACGT\\-]+)$",
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
  for (mut_pattern in names(regex_list)) {
    # Iterate through each (mutually exclusive) mutation regex string
    rgx <- regex_list[[mut_pattern]]
    # Try and match the mutation string to the current regex string
    m <- str_match(mutation, rgx)
    if (all(is.na(m))) {
      # Move to next regex string if no match
      next
    }
    ret$mutation_pattern <- mut_pattern
    rgx_groups <- regex_groups[[mut_pattern]]
    for (field in names(rgx_groups)) {
      # Iterate through the regex capture groups defined for the current mutation pattern
      idx <- rgx_groups[[field]]
      if (any(is.na(idx)) || !is.numeric(idx)) {
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
  # Ensure integer fields (e.g. start and end positions) are stored as integers and not strings
  ret$start <- as.integer(ret$start)
  ret$end <- as.integer(ret$end)
  ret$exon <- as.integer(ret$exon)
  ret$intron_offset_start <- as.integer(ret$intron_offset_start)
  ret$intron_offset_end <- as.integer(ret$intron_offset_end)
  ret$fs_term <- as.integer(ret$fs_term)
  return(ret)
}

parse_aa_change <- function(aa_change, exonic_func) {
  info <- strsplit(aa_change, ":", fixed = TRUE)[[1]]
  info_list <- list(mutation_class = NA, mutation_coding = NA, mutation_protein = NA, nonsynonymous_type = NA)
  info_list$mutation_coding <- grep("^c\\.", info, perl = TRUE, value = TRUE)
  info_list$mutation_protein <- grep("^p\\.", info, perl = TRUE, value = TRUE)
  exon <- gsub("^exon(\\d+)$", "\\1", grep("^exon", info, perl = TRUE, value = TRUE), perl = TRUE)
  if (length(exon) != 1) {
    exon <- NA
  }
  info_list$mutation_coding_info <- parse_mutation(info_list$mutation_coding, aa_change_regex_list$coding, aa_change_regex_groups$coding, exon = exon)
  info_list$mutation_protein_info <- parse_mutation(info_list$mutation_protein, aa_change_regex_list$protein, aa_change_regex_groups$protein, exon = exon)
  if (is.na(info_list$mutation_protein_info$mutation_pattern) ||
      grepl("^(unknown|synonymous.*)$", exonic_func, perl = TRUE)) {
    # No valid AA change, nothing more to do
    return(info_list)
  } else if (grepl("^frameshift", exonic_func, perl = TRUE)) {
    # Framshift flagged by annovar
    info_list$mutation_class <- "frameshift"
  } else if (grepl("^stopgain$", exonic_func, perl = TRUE)) {
    # Stop gain flagged by annovar
    info_list$mutation_class <- "stop_gain"
    if (info_list$mutation_protein_info$mutation_pattern == "sub") {
      # Add additional nonsynonymous type flag for stop gain substitutions
      info_list$nonsynonymous_type <- "sub"
    }
  } else if (is.na(info_list$mutation_protein_info$ref_start) ||
             is.na(info_list$mutation_protein_info$alt) ||
             info_list$mutation_protein_info$ref_start != info_list$mutation_protein_info$alt
             ) {
    # Remaining variants (where ref and alt AAs don't match) are classed as 'nonsynonymous'
    info_list$mutation_class <- "nonsynonymous"
    # Add additional nonsynonymous type flag (sub, dup, del, ins)
    if (info_list$mutation_protein_info$mutation_pattern == "sub") {
      info_list$nonsynonymous_type <- "sub"
    } else if (!is.na(info_list$mutation_protein_info$dupdelins)) {
      info_list$nonsynonymous_type <- info_list$mutation_protein_info$dupdelins
    }
  }
  return(info_list)
}

parse_gene_detail <- function(gene_detail) {
  info <- strsplit(gene_detail, ":", fixed = TRUE)[[1]]
  info_list <- list(mutation = NA, mutation_class = NA)
  info_list$mutation <- grep("^c\\.", info, perl = TRUE, value = TRUE)
  exon <- gsub("^exon(\\d+)$", "\\1", grep("^exon", info, perl = TRUE, value = TRUE), perl = TRUE)
  if (length(exon) != 1) {
    exon <- NA
  }
  info_list$mutation_coding_info <- parse_mutation(info_list$mutation, gene_detail_regex_list, gene_detail_regex_groups, exon = exon)
  if (!is.na(info_list$mutation_coding_info$mutation_pattern)) {
    info_list$mutation_class <- "splice_site"
  }
  return(info_list)
}

parse_chip_def <- function(chip_def, mut_class = NA) {
  info_list <- list()
  info_list$mutation_class <- mut_class
  info_list$mutation <- chip_def
  info_list$mutation_pattern <- NA
  info_list$mutation_info <- parse_mutation(info_list$mutation, chip_def_regex_list, chip_def_regex_groups)
  info_list$mutation_pattern <- info_list$mutation_info$mutation_pattern
  return(info_list)
}

chip_def_match_funcs <- list(
  all = function(...) {
    # Handle "all" CHIP mutation definitions
    args <- list(...)
    
    # Any frameshift, stop_gain, or splice_site c-terminal mutations should FAIL
    if (args$mut_in_c_term &&
        (args$mut_class %in% c("frameshift", "stop_gain", "splice_site"))) {
      return("FAIL")
    }
    
    # Any other variant passes
    return("PASS")
  },
  c_terminal = function(...) {
    # Handle C-terminal mutations
    args <- list(...)
    
    if (args$mut_in_c_term) {
      # C-terminal mutations pass if the CHIP mutation is defined as a C-terminal mutation
      return("PASS")
    } else {
      # Non-C-terminal mutations fail if the CHIP mutation is defined as a C-terminal mutation
      return("FAIL")
    }
  },
  exon = function(...) {
    # Handle CHIP mutations defined by any mutation in a given exon
    args <- list(...)
    chip_info <- args$chip_info
    aa_change_info <- args$aa_change_info
    gene_detail_info <- args$gene_detail_info
    mut_class <- args$mut_class
    
    if (mut_class == "splice_site" && !is.na(gene_detail_info$mutation_coding_info$exon && gene_detail_info$mutation_coding_info$exon == chip_info$mutation_info$exon)) {
      return("PASS")
    } else if (!is.na(aa_change_info$mutation_protein_info$exon) && aa_change_info$mutation_protein_info$exon == chip_info$mutation_info$exon) {
      return("PASS")
    } else {
      return("FAIL")
    }
  },
  exon_aa_insertion = function(...) {
    # Handle CHIP mutations defined by specific amino acid insertions in a given exon
    args <- list(...)
    chip_info <- args$chip_info
    aa_change_info <- args$aa_change_info
    gene_detail_info <- args$gene_detail_info
    mut_class <- args$mut_class
    
    if (!is.na(aa_change_info$mutation_protein_info$exon) && aa_change_info$mutation_protein_info$exon == chip_info$mutation_info$exon) {
      # Mutation has valid exon recorded and matches the CHIP definition exon
      if (!is.na(aa_change_info$mutation_protein_info$dupdelins) && aa_change_info$mutation_protein_info$dupdelins == "ins" && !is.na(aa_change_info$mutation_protein_info$alt) && aa_change_info$mutation_protein_info$alt == chip_info$mutation_info$alt) {
        # Mutation is an insertion and the alternate allele matches the CHIP definition
        return("PASS")
      }
    }
    return("FAIL")
  },
  range = function(...) {
    # Handle amino acid ranges
    args <- list(...)
    chip_info <- args$chip_info
    aa_change_info <- args$aa_change_info
    gene_detail_info <- args$gene_detail_info
    mut_class <- args$mut_class
    
    if (!is.na(aa_change_info$mutation_protein_info$start) &&
        !is.na(aa_change_info$mutation_protein_info$end) &&
        (
          (chip_info$mutation_info$start <= aa_change_info$mutation_protein_info$start && aa_change_info$mutation_protein_info$start <= chip_info$mutation_info$end) ||
          (chip_info$mutation_info$start <= aa_change_info$mutation_protein_info$end && aa_change_info$mutation_protein_info$end <= chip_info$mutation_info$end) ||
          (aa_change_info$mutation_protein_info$start <= chip_info$mutation_info$start && aa_change_info$mutation_protein_info$end >= chip_info$mutation_info$end)
        )) {
      # Passes if mutation amino acid range overlaps with defined CHIP mutation amino acid range
      return("PASS")
    }
    return("FAIL")
  },
  substitution = function(...) {
    # Handle simple substitutions
    args <- list(...)
    chip_info <- args$chip_info
    aa_change_info <- args$aa_change_info
    gene_detail_info <- args$gene_detail_info
    mut_class <- args$mut_class
    
    if (!is.na(aa_change_info$nonsynonymous_type) && aa_change_info$nonsynonymous_type == "sub") {
      if (chip_info$mutation_info$start != aa_change_info$mutation_protein_info$start) {
        return("FAIL")
      }
      if (chip_info$mutation_info$ref_start != aa_change_info$mutation_protein_info$ref_start) {
        return("FAIL")
      }
      if (is.na(chip_info$mutation_info$alt)) {
        return("PASS")
      } else if (chip_info$mutation_info$alt == aa_change_info$mutation_protein_info$alt) {
        return("PASS")
      } else {
        return("FAIL")
      }
    }
    return("FAIL")
  },
  deletion = function(...) {
    # Handle deletions
    args <- list(...)
    chip_info <- args$chip_info
    aa_change_info <- args$aa_change_info
    gene_detail_info <- args$gene_detail_info
    mut_class <- args$mut_class
    
    if (!is.na(aa_change_info$nonsynonymous_type) && aa_change_info$nonsynonymous_type == "del") {
      if (chip_info$mutation_info$start != aa_change_info$mutation_protein_info$start) {
        return("FAIL")
      }
      if (chip_info$mutation_info$ref_start != aa_change_info$mutation_protein_info$ref_start) {
        return("FAIL")
      }
      if (chip_info$mutation_info$end != aa_change_info$mutation_protein_info$end) {
        return("FAIL")
      }
      if (!is.na(aa_change_info$mutation_protein_info$ref_end) && chip_info$mutation_info$ref_end != aa_change_info$mutation_protein_info$ref_end) {
        return("FAIL")
      }
      return("PASS")
    }
    return("FAIL")
  },
  insertion = function(...) {
    # Handle insertions
    args <- list(...)
    chip_info <- args$chip_info
    aa_change_info <- args$aa_change_info
    gene_detail_info <- args$gene_detail_info
    mut_class <- args$mut_class
    
    if (!is.na(aa_change_info$nonsynonymous_type) && aa_change_info$nonsynonymous_type == "ins") {
      if (chip_info$mutation_info$start != aa_change_info$mutation_protein_info$start) {
        return("FAIL")
      }
      if (chip_info$mutation_info$ref_start != aa_change_info$mutation_protein_info$ref_start) {
        return("FAIL")
      }
      if (chip_info$mutation_info$end != aa_change_info$mutation_protein_info$end) {
        return("FAIL")
      }
      if (!is.na(aa_change_info$mutation_protein_info$ref_end) && chip_info$mutation_info$ref_end != aa_change_info$mutation_protein_info$ref_end) {
        return("FAIL")
      }
      if (is.na(chip_info$mutation_info$alt)) {
        return("PASS")
      } else if (chip_info$mutation_info$alt == aa_change_info$mutation_protein_info$alt) {
        return("PASS")
      } else {
        return("FAIL")
      }
    }
    return("FAIL")
  },
  deletion_insertion = function(...) {
    # Handle deletion-insertion events
    args <- list(...)
    chip_info <- args$chip_info
    aa_change_info <- args$aa_change_info
    gene_detail_info <- args$gene_detail_info
    mut_class <- args$mut_class
    
    if (!is.na(aa_change_info$nonsynonymous_type) && aa_change_info$nonsynonymous_type == "delins") {
      if (chip_info$mutation_info$start != aa_change_info$mutation_protein_info$start) {
        return("FAIL")
      }
      if (chip_info$mutation_info$ref_start != aa_change_info$mutation_protein_info$ref_start) {
        return("FAIL")
      }
      if (chip_info$mutation_info$end != aa_change_info$mutation_protein_info$end) {
        return("FAIL")
      }
      if (!is.na(aa_change_info$mutation_protein_info$ref_end) && chip_info$mutation_info$ref_end != aa_change_info$mutation_protein_info$ref_end) {
        return("FAIL")
      }
      if (is.na(chip_info$mutation_info$alt)) {
        return("PASS")
      } else if (chip_info$mutation_info$alt == aa_change_info$mutation_protein_info$alt) {
        return("PASS")
      } else {
        return("FAIL")
      }
    }
    return("FAIL")
  },
  coding_range = function(...) {
    # Handle coding sequence ranges
    args <- list(...)
    chip_info <- args$chip_info
    aa_change_info <- args$aa_change_info
    gene_detail_info <- args$gene_detail_info
    mut_class <- args$mut_class
    
    mut_info <- NA
    if (mut_class == "splice_site" && !all(is.na(gene_detail_info$mutation_coding_info))) {
      mut_info <- gene_detail_info$mutation_coding_info
    } else if (!all(is.na(aa_change_info$mutation_coding_info))) {
      mut_info <- aa_change_info$mutation_coding_info
    } else {
      return("FAIL")
    }
    
    if (is.na(mut_info$start)) {
      return("FAIL")
    }
    if (is.na(mut_info$end)) {
      if (chip_info$mutation_info$start <= mut_info$start && mut_info$start <= chip_info$mutation_info$end) {
        return("PASS")
      } else {
        return("FAIL")
      }
    } else {
      if ((chip_info$mutation_info$start <= mut_info$start && mut_info$start <= chip_info$mutation_info$end) ||
          (chip_info$mutation_info$start <= mut_info$end && mut_info$end <= chip_info$mutation_info$end) ||
          (mut_info$start <= chip_info$mutation_info$start && mut_info$end >= chip_info$mutation_info$end)) {
        return("PASS")
      } else {
        return("FAIL")
      }
    }
    return("FAIL")
  },
  coding_insertion = function(...) {
    # Handle coding insertions
    args <- list(...)
    chip_info <- args$chip_info
    aa_change_info <- args$aa_change_info
    gene_detail_info <- args$gene_detail_info
    mut_class <- args$mut_class
    
    mut_info <- NA
    if (mut_class == "splice_site" && !is.na(gene_detail_info$mutation_coding_info$dupdelins) && gene_detail_info$mutation_coding_info$dupdelins == "ins") {
      mut_info <- gene_detail_info$mutation_coding_info
    } else if (!is.na(aa_change_info$mutation_coding_info$dupdelins) && aa_change_info$mutation_coding_info$dupdelins == "ins") {
      mut_info <- aa_change_info$mutation_coding_info
    } else {
      return("FAIL")
    }
    
    if (is.na(mut_info)) {
      return("FAIL")
    }
    
    if (chip_info$mutation_info$start == mut_info$start && chip_info$mutation_info$end == mut_info$end) {
      return("PASS")
    } else {
      return("FAIL")
    }
    return("FAIL")
  },
  frameshift = function(...) {
    # Handle frameshift events
    args <- list(...)
    chip_info <- args$chip_info
    aa_change_info <- args$aa_change_info
    gene_detail_info <- args$gene_detail_info
    mut_class <- args$mut_class
    
    if (mut_class == "frameshift") {
      if (chip_info$mutation_info$start != aa_change_info$mutation_protein_info$start) {
        return("FAIL")
      }
      if (chip_info$mutation_info$ref_start != aa_change_info$mutation_protein_info$ref_start) {
        return("FAIL")
      }
      return("PASS")
    }
    return("FAIL")
  },
  stop_gain_substitution = function(...) {
    # Handle stop gain events
    args <- list(...)
    chip_info <- args$chip_info
    aa_change_info <- args$aa_change_info
    gene_detail_info <- args$gene_detail_info
    mut_class <- args$mut_class
    
    if (mut_class == "stop_gain" && !is.na(aa_change_info$nonsynonymous_type) && aa_change_info$nonsynonymous_type == "sub") {
      if (chip_info$mutation_info$start != aa_change_info$mutation_protein_info$start) {
        return("FAIL")
      }
      if (chip_info$mutation_info$ref_start != aa_change_info$mutation_protein_info$ref_start) {
        return("FAIL")
      }
      if (is.na(aa_change_info$mutation_protein_info$alt)) {
        return("FAIL")
      }
      if (aa_change_info$mutation_protein_info$alt %in% c("X", "*")) {
        return("PASS")
      } else {
        return("FAIL")
      }
    }
    return("FAIL")
  }
)


# ======================== #
# Main function definition #
# ======================== #

match_mut_def <- function(df, refseq_ensembl_suffix = "ensGene") {
  aachange <- paste("AAChange", refseq_ensembl_suffix, sep = ".")
  genedetail <- paste("GeneDetail", refseq_ensembl_suffix, sep = ".")
  func <- paste("Func", refseq_ensembl_suffix, sep = ".")
  exonic_func <- paste("ExonicFunc", refseq_ensembl_suffix, sep = ".")

  tmp_df <- df
  tmp_df$CHIP_FILTER <- unlist(apply(tmp_df[c(genedetail, aachange, func, exonic_func, "Chr", "Start", "End", "mutation_class", "mutation_definition", "chr", "c_term_genomic_start", "c_term_genomic_end")], 1, function(x) {
    gene_detail <- as.character(x[[1]])
    aa_change <- as.character(x[[2]])
    mut_func <- as.character(x[[3]])
    mut_exonic_func <- as.character(x[[4]])
    mut_chr <- as.character(x[[5]])
    mut_start <- as.numeric(x[[6]])
    mut_end <- as.numeric(x[[7]])
    chip_mut_class <- as.character(x[[8]])
    chip_mut_def <- as.character(x[[9]])
    chip_mut_chr <- as.character(x[[10]])
    chip_mut_c_term_start <- as.numeric(x[[11]])
    chip_mut_c_term_end <- as.numeric(x[[12]])
    
    gene_detail_info <- parse_gene_detail(gene_detail)
    aa_change_info <- parse_aa_change(aa_change, mut_exonic_func)
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
    
    # Rename mutation class for nonsynonymous missense mutations from "nonsynonymous" to "missense"
    # when being compared against a "missense" mutation defined in the CHIP definitions file
    # so it won't be discarded due to a class mismatch in the next step
    if (chip_mut_class == "missense" &&
        mut_class == "nonsynonymous" &&
        !is.na(aa_change_info$nonsynonymous_type) &&
        aa_change_info$nonsynonymous_type == "sub" &&
        !is.na(aa_change_info$mutation_protein_info$ref_start) &&
        !is.na(aa_change_info$mutation_protein_info$alt)) {
      mut_class <- "missense"
    }
    
    # Handle mismatched mutation class
    if (mut_class != chip_mut_class) {
      return("FAIL")
    }
    
    # Determine if mutation is in the c-terminal region
    mut_in_c_term <- FALSE
    if (
      mut_chr == chip_mut_chr && (
        (chip_mut_c_term_start <= mut_start && mut_start <= chip_mut_c_term_end) ||
        (chip_mut_c_term_start <= mut_end && mut_end <= chip_mut_c_term_end) ||
        (mut_start <= chip_mut_c_term_start && mut_end >= chip_mut_c_term_end)
      )
    ) {
      mut_in_c_term <- TRUE
    }
    
    # Call one of the functions stored in the list chip_def_match_funcs.
    # These functions are paired with the mutually-exclusive regex strings defined in chip_def_regex_list,
    # and the name of each regex string and matched function are stored in chip_info$mutation_pattern.
    chip_def_match_mut_pattern_func <- chip_def_match_funcs[[chip_info$mutation_pattern]]
    return(chip_def_match_mut_pattern_func(
      chip_info = chip_info,
      aa_change_info = aa_change_info,
      gene_detail_info = gene_detail_info,
      mut_class = mut_class,
      chip_mut_c_term_start = chip_mut_c_term_start,
      chip_mut_c_term_end = chip_mut_c_term_end,
      mut_start = mut_start,
      mut_end = mut_end,
      mut_in_c_term = mut_in_c_term
    ))
  }))
  
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
