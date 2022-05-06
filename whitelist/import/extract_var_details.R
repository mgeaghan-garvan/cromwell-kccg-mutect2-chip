# =================================================================== #
# Define functions to extract variant information from Annovar output #
# =================================================================== #

# Filter GeneDetail.refGene and AAChange.refGene columns for matching refseq transcript IDs
extractTranscript <- function(column, transcripts, sep) {
  # extractTranscript is to be applied to each row of vars_g_chip_func
  # column: either GeneDetail.refGene or AAChange.refGene entry for a given row of vars_g_chip_func
  # transcripts: a comma-delimited list of representative transcripts to match against
  transcripts_list <- strsplit(transcripts, ",")[[1]]
  column_list <- strsplit(column, sep)[[1]]
  column_match <- list()
  if (length(column_list) == 0) {
    return("nan")
  }
  for(i in 1:length(column_list)) {
    column_match[[i]] <- TRUE %in% (transcripts_list %in% strsplit(column_list, ":")[[i]])
  }
  column_match <- unlist(column_match)
  if(TRUE %in% column_match) {
    return(paste(column_list[column_match], collapse = sep))
  } else {
    return("nan")
  }
}

extractNonsyn <- function(AAChange){
  # extractNonsyn is to be applied to each row of vars_g_chip_func
  # AAChange: the AAChange.transcript entry for a given row of vars_g_chip_func
  AAChange_list <- strsplit(AAChange, ",")[[1]]
  protChange_list <- unlist(lapply(
    strsplit(AAChange_list, ":"), function(x) { grep("p.", x, value = T, fixed = T) }
  ))
  if(length(protChange_list) > 0) {
    return(paste(gsub("p\\.", "", protChange_list), collapse = ","))
  } else {
    return("nan")
  }
}

extractNonsynExon <- function(AAChange){
  # extractNonsyn is to be applied to each row of vars_g_chip_func
  # AAChange: the AAChange.transcript entry for a given row of vars_g_chip_func
  AAChange_list <- strsplit(AAChange, ",")[[1]]
  protChange_list <- unlist(lapply(
    strsplit(AAChange_list, ":"), function(x) { grep("exon", x, value = T, fixed = T) }
  ))
  if(length(protChange_list) > 0) {
    return(paste(gsub("exon", "", protChange_list), collapse = ","))
  } else {
    return("nan")
  }
}

extractSpliceExons <- function(GeneDetail) {
  # extractSpliceExons is to be applied to each row of vars_g_chip_func
  # GeneDetail: the GeneDetail.refGene entry for a given row of vars_g_chip_func
  GeneDetail_list <- strsplit(GeneDetail, ";")[[1]]
  change_list <- unlist(lapply(
    strsplit(GeneDetail_list, ":"), function(x) { grep("exon", x, value = T, fixed = T) }
  ))
  if(length(change_list) > 0) {
    return(paste(gsub("exon", "", change_list), collapse = ";"))
  } else {
    return("nan")
  }
}