# Get sequence context around variants and set filter to manual review if AD < 10 or VAF < 0.1
get_hprs <- function(seq, min_size = 5) {
  # Sliding window algorithm to pick out all homopolymer regions within a sequence
  hpr_list <- list()
  hidx <- 1
  l <- nchar(seq)
  i <- 1
  while (i <= (l - min_size + 1)) {
    b1 <- substr(seq, i, i)
    j <- i + 1
    while (j <= l) {
      b2 <- substr(seq, j, j)
      if (b2 == b1) {
        j <- j + 1
      } else {
        break
      }
    }
    subseq_len <- j - i
    if (subseq_len >= min_size) {
      hpr_list[[hidx]] <- data.frame(start = i, end = (j - 1), length = (j - i), base = b1)
      hidx <- hidx + 1
    }
    i <- j
  }
  return(do.call(rbind, hpr_list))
}

apply_homopolymer_filter <- function(df, fasta_file, HP_SNVS = FALSE) {
  # Should we look at homopolymer SNVs as well as INDELs?
  # NOTE: For now, we'll just be doing INDELs (HP_SNVS = FALSE), but in the future we should make this configurable via CLI arguments
  var_selection <- NA
  if (HP_SNVS) {
    # Get all variants passing previous filters
    var_selection <- df$COMBINED_FILTER == "PASS"
  } else {
    # Just retrieve INDELs
    var_selection <- grepl("(insertion|deletion)", df$ExonicFunc.refGene, perl = TRUE)
  }
  
  # Create sub-set of variants that pass previous filters. Get basic variant information
  vars_pass <- df[var_selection, c("Chr", "POS", "REF", "ALT")]
  # Construct bed-formatted variant data-frame
  vars_pass$Start <- vars_pass$POS
  vars_pass$End <- vars_pass$Start + nchar(vars_pass$REF) - 1
  vars_pass_bed <- unique(vars_pass[, c("Chr", "Start", "End")])
  vars_pass_bed$Start = as.integer(vars_pass_bed$Start) - 11  # Zero-based coordinates for BED format, get 10bp upstream
  vars_pass_bed$End = as.integer(vars_pass_bed$End) + 10  # Get 10bp downstream
  # Write to a bed file and run bedtools getfasta to get sequence context, then read back in
  write.table(vars_pass_bed, "_tmp_pass.bed", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
  system(paste("bedtools getfasta -fi ", fasta_file, " -bed _tmp_pass.bed -fo _tmp_pass.fa", sep = ""))
  vars_pass_fa <- scan("_tmp_pass.fa", character())
  if (!(length(vars_pass_fa) %% 2 == 0)) {
    # Something went wrong and the fasta file is improperly formatted (number of lines not a mulitple of 2)
    stop("ERROR: Could not successfully retrieve variant sequence contexts.")
  } else {
    pass_seq = list()
    i = 1
    j = 1
    for (s in vars_pass_fa) {
      if (i %% 2 == 1 && !grepl("^>", s, perl = TRUE)) {
        # Second of every pair of lines should not start with a ">" (should be the sequence itself)
        stop("ERROR: Could not successfully retrieve variant sequence contexts.")
      } else if (i %% 2 == 0 && grepl("^>", s, perl = TRUE)) {
        # First of every pair of lines should start with ">" (this is the sequence name)
        stop("ERROR: Could not successfully retrieve variant sequence contexts.")
      } else {
        if (i %% 2 == 1) {
          # Parse the sequence name: retrieve chr, start, and end positions
          name <- s
          chr <- gsub("\\:.*", "", name, perl = TRUE)
          chr <- gsub("^>", "", chr, perl = TRUE)
          start <- gsub(".*\\:", "", name, perl = TRUE)
          start <- gsub("\\-.*", "", start, perl = TRUE)
          start <- as.integer(start) + 11
          end <- gsub(".*\\-", "", name, perl = TRUE)
          end <- as.integer(end) - 10
          pass_seq[[j]] <- c(chr, start, end, "")
          i <- i + 1
        } else if (i %% 2 == 0) {
          # Get the sequence context
          pass_seq[[j]][4] = s
          i <- i + 1
          j <- j + 1
        }
      }
    }
    # Merge results into a single dataframe, rename the columns, and format the position columns as integers
    pass_seq <- as.data.frame(do.call(rbind, pass_seq))
    colnames(pass_seq) <- c("Chr", "Start", "End", "VAR_SEQ_CONTEXT")
    vars_pass$Start <- as.integer(vars_pass$Start)
    vars_pass$End <- as.integer(vars_pass$End)
    vars_pass <- merge(vars_pass, pass_seq, all.x = TRUE)
    # Find homopolymer regions
    vars_pass_filter <- do.call(rbind, apply(vars_pass[c("Start", "End", "REF", "ALT", "VAR_SEQ_CONTEXT")], 1, function(x) {
      # Get alternate allele sequence context
      s <- as.integer(x[[1]])
      e <- as.integer(x[[2]])
      r <- as.character(x[[3]])
      a <- strsplit(as.character(x[[4]]), ",")[[1]]
      l <- e - s + 1  # l <- nchar(r)
      rs <- as.character(x[[5]])
      if (length(a) == 0 || is.null(a) || is.na(a)) {
        return(data.frame(VAR_ALT_SEQ_CONTEXT = "", HOMOPOLYMER_FILTER = "FAIL"))  # This shouldn't ever occur, but if it does we have an invalid variant so it should FAIL
      }
      as_list = list()
      filter_list = list()
      for (allele_idx in 1:length(a)) {
        a_i <- a[allele_idx]
        as <- paste(
          substr(rs, 1, 10),
          a_i,
          substr(rs, (11 + l), nchar(rs)),
          sep = ""
        )
        # Find all homopolymers in both ref and alt sequences
        hpr_ref <- get_hprs(rs)
        hpr_alt <- get_hprs(as)
        filter <- ""
        # Find variants in homopolymer regions
        if (is.null(hpr_ref) && is.null(hpr_alt)) {
          # No homopolymers - filter passes
          filter <- "PASS"
        } else if (is.null(hpr_ref) && !is.null(hpr_alt)) {
          # Alternate allele creates homopolymer
          filter <- "HOMOPOLYMER_VARIANT"
        } else if (!is.null(hpr_ref) && is.null(hpr_alt)) {
          # Alternate allele destroys homopolymer
          filter <- "HOMOPOLYMER_VARIANT"
        } else if (!identical(hpr_ref[c("length", "base")], hpr_alt[c("length", "base")])) {
          # Reference and alternate sequences have differing homopolymers (i.e. variant creates/destroys a homopolymer sequence)
          filter <- "HOMOPOLYMER_VARIANT"
        } else if (all(hpr_ref$end <= 10) && all(hpr_alt$end <= 10)) {
          # All homopolymers reside within first 10 bp (outside of variant range which starts at bp 11)
          filter <- "PASS"
        } else if (l == 1) {
          # All remaining insertions should pass
          filter <- "PASS"
        } else {
          # All remaining variants are deletions
          filter <- "PASS"
          s_del <- 11
          e_del <- 11 + l - 1
          for (i in 1:dim(hpr_ref)[1]) {
            s_hpr <- hpr_ref[i, "start"]
            e_hpr <- hpr_ref[i, "end"]
            if (
              ((s_hpr <= s_del) && (s_del <= e_hpr)) ||
              ((s_hpr <= e_del) && (e_del <= e_hpr)) ||
              ((s_del < s_hpr) && (e_del > e_hpr))
            ) {
              # Deletion starts or ends within, or contains a homopolymer region
              filter <- "HOMOPOLYMER_VARIANT"
              break
            }
          }
        }
        as_list[[allele_idx]] <- as
        filter_list[[allele_idx]] <- filter
      }
      return(data.frame(VAR_ALT_SEQ_CONTEXT = paste(as_list, collapse = ","), HOMOPOLYMER_FILTER = paste(filter_list, collapse = ",")))
    }))
    vars_pass$VAR_ALT_SEQ_CONTEXT <- vars_pass_filter$VAR_ALT_SEQ_CONTEXT
    vars_pass$HOMOPOLYMER_FILTER <- vars_pass_filter$HOMOPOLYMER_FILTER
  }
  
  # Discard the temporary 'Start' and 'End' columns
  vars_pass <- vars_pass[!(colnames(vars_pass) %in% c("Start", "End"))]
  preserve_col_order <- union(colnames(df), colnames(vars_pass))
  # Merge the new sequence and filter columns into the main dataframe
  df <- merge(df, vars_pass, all.x = TRUE)
  df <- df[preserve_col_order]
  # Apply homopolymer filters
  df$HOMOPOLYMER_FILTER[is.na(df$HOMOPOLYMER_FILTER)] <- "PASS"
  df$HOMOPOLYMER_FILTER <- unlist(apply(df[c("HOMOPOLYMER_FILTER", "AD", "VAF")], 1, function(x) {
    filter <- as.character(x[[1]])
    if (filter == "PASS") {
      return("PASS")
    }
    # Get filter, AD, and VAF values
    # Multiallelic variants will have multiple comma-separated values for each
    filter <- strsplit(filter, ",")[[1]]
    ad <- as.integer(strsplit(as.character(x[[2]]), ",")[[1]])
    ad <- ad[2:length(ad)]
    vaf <- as.numeric(strsplit(as.character(x[[3]]), ",")[[1]])
    if (length(filter) != length(ad) && length(filter) != length(vaf)) {
      # Something's not right with the filter, AD, or VAF columns - investigate
      # Hopefully this never happens!
      return("MANUAL_REVIEW")
    }
    # For each alternate allele, check AD and VAF against filter thresholds
    new_filter_list <- list()
    for (filter_idx in 1:length(filter)) {
      filter_i <- filter[filter_idx]
      ad_i <- ad[filter_idx]
      vaf_i <- vaf[filter_idx]
      if (filter_i == "HOMOPOLYMER_VARIANT" & (ad_i < 10 || vaf_i < 0.1)) {
        new_filter_list[[filter_idx]] <- "FAIL"
      } else {
        new_filter_list[[filter_idx]] <- filter_i
      }
    }
    # Update filter column
    if (all(new_filter_list %in% c("PASS", "HOMOPOLYMER_VARIANT"))) {
      return("PASS")
    } else if (all(new_filter_list == "FAIL")) {
      return("FAIL")
    } else {
      return("MANUAL_REVIEW")
    }
  }))
  
  # Update combined filter column
  df$COMBINED_FILTER <- apply(df[c("COMBINED_FILTER", "HOMOPOLYMER_FILTER")], 1, function(x) {
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
  return(unique(df))
}