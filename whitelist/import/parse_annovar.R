parse_annovar <- function(df, vars_variant_func, ensGene = TRUE, refGene = FALSE, filter = c("exonic", "splicing")) {
  col_suffix <- ""
  trans_regex <- ""
  if (ensGene) {
    col_suffix <- "ensGene"
    trans_regex <- "^(ENST\\d+)(\\.\\d+)?$"
  } else if (refGene) {
    col_suffix <- "refGene"
    trans_regex <- "^(N[MR]_\\d+)(\\.\\d+)?$"
  } else {
    warning("Neither ensGene or refGene annotations selected, returning original dataframe")
    return(df)
  }
  
  aachange <- paste("AAChange", col_suffix, sep = ".")
  genedetail <- paste("GeneDetail", col_suffix, sep = ".")
  func <- paste("Func", col_suffix, sep = ".")
  transcript <- paste("Transcript", col_suffix, sep = ".")
  
  # Split rows by RefSeq or Ensembl IDs in both GeneDetail (noncoding variant annotations, e.g. splicing variants) and AAChange (coding variant annotations)
  if (ensGene) {
    vars_g <- as.data.frame(separate_rows(df, GeneDetail.ensGene, sep = ";"))
    vars_g <- as.data.frame(separate_rows(vars_g, AAChange.ensGene, sep = ","))
    vars_g <- as.data.frame(separate_rows(vars_g, AAChange.ensGene, sep = ";"))
  } else {
    vars_g <- as.data.frame(separate_rows(df, GeneDetail.refGene, sep = ";"))
    vars_g <- as.data.frame(separate_rows(vars_g, AAChange.refGene, sep = ","))
    vars_g <- as.data.frame(separate_rows(vars_g, AAChange.refGene, sep = ";"))
  }
  
  # Filter variants for exonic or splicing variants
  if (!is.na(filter) && !is.null(filter)) {
    vars_g_filter <- vars_g[grepl(paste("(", paste(filter, collapse = "|"), ")", sep = ""), vars_g[[func]], perl = TRUE),]
  } else {
    vars_g_filter <- vars_g
  }
  
  # Drop potential duplicates
  vars_g_filter <- unique(vars_g_filter)
  
  # Construct a column with the RefSeq or Ensembl IDs from both GeneDetail and AAChange columns, and split rows again on this column
  vars_g_func_acc <- unlist(apply(vars_g_filter[c(aachange, genedetail)], 1, function(x) {
    aa <- as.character(x[[1]])
    gd <- as.character(x[[2]])
    aagd <- paste(aa, gd, sep = ":")
    ann <- strsplit(aagd, ":")[[1]]
    acc <- grep(trans_regex, ann, perl = TRUE, value = TRUE)
    acc <- gsub(trans_regex, "\\1", acc, perl = TRUE)
    return(paste(acc, collapse = ";"))
  }))
  vars_g_filter[[transcript]] <- vars_g_func_acc
  
  if (ensGene) {
    vars_g_filter <- as.data.frame(separate_rows(vars_g_filter, Transcript.ensGene, sep = ";"))
  } else {
    vars_g_filter <- as.data.frame(separate_rows(vars_g_filter, Transcript.refGene, sep = ";"))
  }
  
  # Define function to return matching transcript information from AAChange/GeneDetail column
  acc_func <- function(x) {
    acc <- as.character(x[[1]])
    ann <- as.character(x[[2]])
    ret <- grep(paste("(^|:)", acc, "(\\.\\d+)?(:|$)", sep = ""), ann, perl = TRUE, value = TRUE)
    if (length(ret) == 0) {
      return("")
    } else if (length(ret) > 1) {
      warning(paste("Multiple accession numbers found in record:", ann, "- skipping..."))
      return("")
    } else {
      return(ret)
    }
  }
  
  vars_g_filter[[aachange]] <- unlist(apply(vars_g_filter[c(transcript, aachange)], 1, acc_func))
  vars_g_filter[[genedetail]] <- unlist(apply(vars_g_filter[c(transcript, genedetail)], 1, acc_func))
  
  # Drop potential duplicates and rows without a RefSeq ID
  vars_g_filter <- unique(vars_g_filter[vars_g_filter[[transcript]] != "",])
  
  # Split rows once again on the Func.refGene/Func.ensGene column
  if (ensGene) {
    vars_g_filter <- as.data.frame(separate_rows(vars_g_filter, Func.ensGene, sep = ";"))
  } else {
    vars_g_filter <- as.data.frame(separate_rows(vars_g_filter, Func.refGene, sep = ";"))
  }
  
  # Merge with the variant_func dataframe
  vars_g_filter <- merge(vars_g_filter, vars_variant_func, by.x = c(func, transcript, "Chr", "Start", "End", "Ref", "Alt"), by.y = c("Func", "Transcript", "Chr", "Start", "End", "Ref", "Alt"), all.x = TRUE)
  
  # Filter variants again for exonic or splicing variants
  if (!is.na(filter) && !is.null(filter)) {
    vars_g_filter <- vars_g_filter[grepl(paste("^(", paste(filter, collapse = "|"), ")$", sep = ""), vars_g_filter[[func]], perl = TRUE),]
  }
  
  # Drop potential duplicates again
  vars_g_filter <- unique(vars_g_filter)
  
  return(vars_g_filter)
}