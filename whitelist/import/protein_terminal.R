get_protein_terminal_status <- function(df, transcript_prot_file, parse_aa_change_func) {
  prot_lengths <- read.csv(transcript_prot_file, header = TRUE)
  colnames(prot_lengths) <- c("refseq_mrna", "hgnc_symbol", "prot_length")
  df$AAChange.N_TERM_10PCT <- unlist(lapply(df$AAChange.refGene, function(x) {
    t_list <- strsplit(x, ",")[[1]]
    t_list <- strsplit(t_list, ":")
    ret <- list()
    if (length(t_list) == 0) {
      return("")
    }
    for (i in 1:length(t_list)) {
      t <- t_list[[i]]
      t_in_prot_lengths <- t %in% prot_lengths$refseq_mrna
      if (!any(t_in_prot_lengths)) {
        ret[[i]] <- "NA"
      } else {
        t_name <- t[t_in_prot_lengths][1]
        l <- prot_lengths$prot_length[prot_lengths$refseq_mrna == t_name][1]
        a <- grep("^p\\.", t, perl = TRUE, value = TRUE)
        if (length(a) == 0) {
          ret[[i]] <- "NA"
          next
        }
        a <- gsub("^p\\.", "", a, perl = TRUE)
        p <- parse_aa_change_func(a)$start_pos
        if (is.na(p)) {
          ret[[i]] <- "NA"
        } else if (p > l) {
          ret[[i]] <- NA
        } else {
          ret[[i]] <- (p/l) < 0.1
        }
      }
    }
    ret <- as.character(unlist(ret))
    return(paste(ret, collapse = ","))
  }))
  df$AAChange.C_TERM_10PCT <- unlist(lapply(df$AAChange.refGene, function(x) {
    t_list <- strsplit(x, ",")[[1]]
    t_list <- strsplit(t_list, ":")
    ret <- list()
    if (length(t_list) == 0) {
      return("")
    }
    for (i in 1:length(t_list)) {
      t <- t_list[[i]]
      t_in_prot_lengths <- t %in% prot_lengths$refseq_mrna
      if (!any(t_in_prot_lengths)) {
        ret[[i]] <- "NA"
      } else {
        t_name <- t[t_in_prot_lengths][1]
        l <- prot_lengths$prot_length[prot_lengths$refseq_mrna == t_name][1]
        a <- grep("^p\\.", t, perl = TRUE, value = TRUE)
        if (length(a) == 0) {
          ret[[i]] <- "NA"
          next
        }
        a <- gsub("^p\\.", "", a, perl = TRUE)
        p <- parse_aa_change_func(a)$start_pos
        if (is.na(p)) {
          ret[[i]] <- "NA"
        } else if (p > l) {
          ret[[i]] <- NA
        } else {
          ret[[i]] <- (p/l) > 0.9
        }
      }
    }
    ret <- as.character(unlist(ret))
    return(paste(ret, collapse = ","))
  }))
}