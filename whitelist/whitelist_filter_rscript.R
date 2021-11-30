# library(data.table, quietly = T)
library(tidyr)
library(stringr)

# Check command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop(paste(
    "Incorrect number of arguments!",
    "Usage: Rscript whitelist_filter_rscript.R <ANNOVAR OUTPUT TABLE> <ANNOVAR OUTPUT VCF> <TUMOR SAMPLE NAME> <GNOMAD SOURCE> <GNOMAD SUBPOPULATION CODE> <TREAT MISSING AF AS RARE> <CHIP DEFINITION FILE>",
    "    gnomAD source:              'exome', 'genome', or 'both'",
    "    gnomAD subpopulation codes: 'AF' (all), 'AF_afr', 'AF_sas', 'AF_amr', 'AF_eas', 'AF_nfe', 'AF_fin', 'AF_asj'",
    "    TREAT MISSING AF AS RARE:   'TRUE' = variants not annotated in gnomAD are assumed to be rare and are given an allele frequency of 0; 'FALSE' = variants not annotated in gnomAD will not pass the gnomAD hard filter.",
    "    CHIP DEFINITION FILE:       csv file containing CHIP variant definitions",
    sep = "\n"))
}
annovar_text_out <- args[1]
annovar_vcf_out <- args[2]
tumor_sample_name <- args[3]
gnomad_source <- args[4]
gnomad_pop <- args[5]
treat_missing_as_rare <- args[6]
chip_def_file <- args[7]

if (!file.exists(annovar_text_out)) {
  stop("Input annovar table file does not exist.")
}
if (!file.exists(annovar_vcf_out)) {
  stop("Input annovar vcf file does not exist.")
}
if (!file.exists(chip_def_file)) {
  stop("Chip variant definition file does not exist.")
}
if (!(gnomad_source %in% c("exome", "genome", "both"))) {
  stop("Invalid gnomAD source.")
}
if (!(gnomad_pop %in% c("AF", "AF_afr", "AF_sas", "AF_amr", "AF_eas", "AF_nfe", "AF_fin", "AF_asj"))) {
  stop("Invalid gnomAD population.")
}
if (!(treat_missing_as_rare %in% c("TRUE", "FALSE"))) {
  stop("Positional argument #4 'TREAT MISSING AF AS RARE' must be 'TRUE' or 'FALSE'")
} else {
  treat_missing_as_rare <- (treat_missing_as_rare == "TRUE")
}
annovar_text_out_regex <- "_multianno\\.txt$"
if (!grepl(annovar_text_out_regex, annovar_text_out, perl = TRUE)) {
  stop("Invalid input annovar table file. Must be an annovar text output ('*_multianno.txt').")
}
annovar_vcf_out_regex <- "_multianno\\.vcf$"
if (!grepl(annovar_vcf_out_regex, annovar_vcf_out, perl = TRUE)) {
  stop("Invalid input annovar vcf file. Must be an annovar vcf output ('*_multianno.vcf').")
}

# Load CHIP variant/gene lists
chip_vars <- read.csv(chip_def_file)

# Define sample ID
sample_id <- gsub(annovar_text_out_regex, "", annovar_text_out, perl = TRUE)
sample_id <- gsub("^.*\\/([^\\/]+)$", "\\1", sample_id, perl = TRUE)

# Load annovar variant annotations
vars <- read.table(annovar_text_out, sep = "\t", header = TRUE)
vars$Sample <- sample_id

# Load annovar vcf file
vcf <- scan(annovar_vcf_out, character(), sep = "\n")
vcf_header <- grep("^#CHROM", vcf, perl = TRUE, value = TRUE)
vcf_header <- gsub("^#", "", vcf_header, perl = TRUE)
vcf_header <- strsplit(vcf_header, "\t")[[1]]

# Ensure the following columns are present
req_cols <- c("Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene", "GeneDetail.refGene", "ExonicFunc.refGene", "AAChange.refGene")
if (gnomad_source == "both") {
  gnomad_g_cols <- c("AF", "AF_popmax", "AF_male", "AF_female", "AF_raw", "AF_afr", "AF_sas", "AF_amr", "AF_eas", "AF_nfe", "AF_fin", "AF_asj", "AF_oth", "non_topmed_AF_popmax", "non_neuro_AF_popmax", "non_cancer_AF_popmax", "controls_AF_popmax")
  gnomad_e_cols <- c("AF.1", "AF_popmax.1", "AF_male.1", "AF_female.1", "AF_raw.1", "AF_afr.1", "AF_sas.1", "AF_amr.1", "AF_eas.1", "AF_nfe.1", "AF_fin.1", "AF_asj.1", "AF_oth.1", "non_topmed_AF_popmax.1", "non_neuro_AF_popmax.1", "non_cancer_AF_popmax.1", "controls_AF_popmax.1")
  req_cols <- c(req_cols, gnomad_g_cols, gnomad_e_cols)
} else if (gnomad_source == "exome" || gnomad_source == "genome") {
  gnomad_cols <- c("AF", "AF_popmax", "AF_male", "AF_female", "AF_raw", "AF_afr", "AF_sas", "AF_amr", "AF_eas", "AF_nfe", "AF_fin", "AF_asj", "AF_oth", "non_topmed_AF_popmax", "non_neuro_AF_popmax", "non_cancer_AF_popmax", "controls_AF_popmax")
  req_cols <- c(req_cols, gnomad_cols)
} else {
  stop("Invalid gnomAD source.")
}
stopifnot(all(req_cols %in% colnames(vars)))

# Rename Otherinfo columns
# NOTE: These are hard-coded to match the output of annovar as of 2021-11-22
#       When converting from vcf to avinput format, the following command is used by table_annovar.pl:
#       convert2annovar.pl -includeinfo -allsample -withfreq -format vcf4 INPUT.vcf > OUTPUT.avinput
#       This produces 3 "Otherinfo" columns corresponding to total allele frequency across all samples, quality score, and read deapth
#       The following 9+N columns (for N samples) correspond to the input vcf file's columns
otherinfo_cols <- c("ANNOVAR_alt_af", "ANNOVAR_qual", "ANNOVAR_alt_ad", vcf_header)
colnames(vars)[grepl("Otherinfo", colnames(vars), fixed = TRUE)] <- otherinfo_cols

# ===== PRE-FILTERING =====

# Get gnomAD allele frequencies
if (gnomad_source == "both") {
  new_gnomad_g_cols <- paste("gnomAD_genome_", gnomad_g_cols, sep = "")
  new_gnomad_e_cols <- gsub("\\.\\d+$", "", gnomad_e_cols, perl = TRUE)
  new_gnomad_e_cols <- paste("gnomAD_exome_", new_gnomad_e_cols, sep = "")

  colnames(vars)[colnames(vars) %in% gnomad_g_cols] <- new_gnomad_g_cols
  colnames(vars)[colnames(vars) %in% gnomad_e_cols] <- new_gnomad_e_cols

  gnomad_pop_columns <- paste("gnomAD", c("genome", "exome"), gnomad_pop, sep = "_")

  vars$gnomAD_AF <- apply(vars[gnomad_pop_columns], 1, function(x) {
    af_g <- x[[1]]
    af_e <- x[[2]]
    missing_vals <- c(".", "", NA)
    if(!(af_e %in% missing_vals)) {
      return(as.numeric(af_e))
    } else if(!(af_g %in% missing_vals)) {
      return(as.numeric(af_g))
    } else if(treat_missing_as_rare) {
      return(0)
    } else {
      return(NA)
    }
  })
} else {
  new_gnomad_cols <- paste("gnomAD_", gnomad_source, "_", gnomad_cols, sep = "")

  colnames(vars)[colnames(vars) %in% gnomad_cols] <- new_gnomad_cols

  gnomad_pop_columns <- paste("gnomAD", gnomad_source, gnomad_pop, sep = "_")

  vars$gnomAD_AF <- apply(vars[gnomad_pop_columns], 1, function(x) {
    af <- x[[1]]
    missing_vals <- c(".", "", NA)
    if(!(af %in% missing_vals)) {
      return(as.numeric(af))
    } else if(treat_missing_as_rare) {
      return(0)
    } else {
      return(NA)
    }
  })
}


# Filter by AD, DP, AF, F1R2/F2R1 VCF fields and gnomAD frequency
get_format_field <- function(x, format_field) { grep(paste("^", x, "$", sep = ""), format_field, perl = TRUE) }
vars$HARD_FILTER <- apply(vars[c("FORMAT", tumor_sample_name, "gnomAD_AF")], 1, function(x) {
  format_field <- strsplit(x[[1]], ":")[[1]]
  sample_field <- strsplit(x[[2]], ":")[[1]]
  gnomad_af <- as.numeric(x[[3]])
  # Allelic depth
  var_n_ad_i <- get_format_field("AD", format_field)
  var_n_ad <- as.integer(strsplit(sample_field[var_n_ad_i], ",")[[1]])
  var_n_ad_gte_3 <- var_n_ad >= 3
  # Site read depth
  var_n_dp_i <- get_format_field("DP", format_field)
  var_n_dp <- as.integer(sample_field[var_n_dp_i])
  var_n_dp_gte_20 <- var_n_dp >= 20
  # VAF
  var_n_vaf_i <- get_format_field("AF", format_field)
  var_n_vaf <- as.numeric(strsplit(sample_field[var_n_vaf_i], ",")[[1]])
  var_n_vaf_gte_2pc <- var_n_vaf >= 0.02
  var_n_vaf_lt_35pc <- var_n_vaf < 0.35
  # Forward and reverse stand support
  var_n_for_i <- get_format_field("F1R2", format_field)
  var_n_for <- as.integer(strsplit(sample_field[var_n_for_i], ",")[[1]])
  var_n_rev_i <- get_format_field("F2R1", format_field)
  var_n_rev <- as.integer(strsplit(sample_field[var_n_rev_i], ",")[[1]])
  var_n_for_rev_gte_1 <- all(var_n_for >= 1) && all(var_n_rev >= 1)
  # gnomAD frequency
  var_n_gnomad_af_lt_0.1pc <- !is.na(gnomad_af) & (gnomad_af < 0.001)
  filter_pass <- (
    all(var_n_ad_gte_3[2:length(var_n_ad_gte_3)]) &  # first allele = ref
      var_n_dp_gte_20 &
      all(var_n_vaf_gte_2pc) &
      all(var_n_vaf_lt_35pc) &
      var_n_for_rev_gte_1 &
      var_n_gnomad_af_lt_0.1pc
  )
  filter_manual_review <- (
    any(var_n_ad_gte_3[2:length(var_n_ad_gte_3)]) &  # first allele = ref; if any, but not all pass, go to manual review
      var_n_dp_gte_20 &
      any(var_n_vaf_gte_2pc) &  # some but not all are passing
      var_n_for_rev_gte_1 &
      var_n_gnomad_af_lt_0.1pc
  )
  if(filter_pass) {
    return("PASS")
  } else if(filter_manual_review) {
    return("MANUAL_REVIEW")
  } else {
    return("FAIL")
  }
})
vars$COMB_FILTER = apply(vars[c("FILTER", "HARD_FILTER")], 1, function(x) {
  if(x[[1]] == "PASS") {
    return(x[[2]])
  } else {
    return("FAIL")
  }
})

# ===== PARSE ANNOVAR OUTPUT =====

# Split multi-gene variants
vars_g <- as.data.frame(separate_rows(vars, Gene.refGene, sep = ";"))
vars_g$AAChange.refGene <- apply(vars_g[, c("Gene.refGene", "AAChange.refGene")], 1, function(row) {
  a <- strsplit(row[[2]], ",")[[1]]
  aa <- a[grepl(paste("(^|:)", row[[1]], "(:|$)", sep = ""), a, perl = TRUE)]
  aaa <- paste(aa, collapse = ",")
  return(aaa)
})

# Filter for CHIP genes
vars_g_chip <- vars_g[vars_g$Gene.refGene %in% chip_vars$Gene,]

# Filter variants for exonic or splicing variants
vars_g_chip_func <- vars_g_chip[
  (
    grepl("exonic", vars_g_chip$Func.refGene, fixed = T) |
      grepl("splicing", vars_g_chip$Func.refGene, fixed = T)
  ),
]

# Merge filtered variants and CHIP gene annotation data frames
vars_g_chip_func <- merge(vars_g_chip_func, chip_vars, by.x = "Gene.refGene", by.y = "Gene")

# Filter GeneDetail.refGene and AAChange.refGene columns for matching refseq transcript IDs
extractTranscript <- function(column, transcripts, sep) {
  # extractTranscript is to be applied to each row of vars_g_chip_func
  # column: either GeneDetail.refGene or AAChange.refGene entry for a given row of vars_g_chip_func
  # transcripts: a comma-delimited list of representative transcripts to match against
  transcripts_list <- strsplit(transcripts, ",")[[1]]
  column_list <- strsplit(column, sep)[[1]]
  column_match <- list()
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

vars_g_chip_func$GeneDetail.transcript <- unlist(apply(vars_g_chip_func[, c("GeneDetail.refGene", "Transcript_Accession")], 1, function(x) {
  extractTranscript(x[1], x[2], ";")
}))

vars_g_chip_func$AAChange.transcript <- unlist(apply(vars_g_chip_func[, c("AAChange.refGene", "Transcript_Accession")], 1, function(x) {
  extractTranscript(x[1], x[2], ",")
}))

vars_g_chip_func$AAChange.protChange <- unlist(lapply(vars_g_chip_func[, c("AAChange.transcript")], extractNonsyn))

vars_g_chip_func$AAChange.exon <- unlist(lapply(vars_g_chip_func[, c("AAChange.transcript")], extractNonsynExon))

vars_g_chip_func$GeneDetail.exon <- unlist(lapply(vars_g_chip_func[, c("GeneDetail.transcript")], extractSpliceExons ))

# Remove variants with no matching affected transcript
vars_g_chip_func_filtered <- vars_g_chip_func[!(
  vars_g_chip_func$GeneDetail.transcript == "nan" &
    vars_g_chip_func$AAChange.transcript == "nan" &
    vars_g_chip_func$AAChange.protChange == "nan" &
    vars_g_chip_func$AAChange.exon == "nan" &
    vars_g_chip_func$GeneDetail.exon == "nan" &
    !(vars_g_chip_func$Nonsynonymous %in% c("any", "c_term")) &
    !(vars_g_chip_func$Frameshift %in% c("any", "c_term")) &
    !(vars_g_chip_func$Stop_Gain %in% c("any", "c_term")) &
    !(vars_g_chip_func$Splicing %in% c("any", "c_term")) &
    !(vars_g_chip_func$Putative_Nonsynonymous %in% c("any", "c_term")) &
    !(vars_g_chip_func$Putative_Frameshift %in% c("any", "c_term")) &
    !(vars_g_chip_func$Putative_Stop_Gain %in% c("any", "c_term")) &
    !(vars_g_chip_func$Putative_Splicing %in% c("any", "c_term"))
), ]

# Functions for parsing condition and AA change codes
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

# Functions for matching variants with CHIP conditions
match_ns_conditions <- function(conditions, aa_changes, aa_exons) {
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
    }
    cond_details <- parse_cond(cond)
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
      a_details <- parse_aa_change(a)
      if(all(is.na(a_details))) {
        next
      }
      if(a_details$mut_type == "aa_change" & (a_details$start_pos %in% cond_details$start_pos:cond_details$end_pos)) {
        # Current variant is a single AA substitution and start position is in defined range
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
      } else if(a_details$mut_type == "ins" & (a_details$start_pos %in% cond_details$start_pos:cond_details$end_pos)) {
        # Current variant is an insertion and start position is in defined range
        if(cond_details$cond_type == "aa_range") {
          # Current variant in defined range
          whitelist <- TRUE
          break
        } else if(cond_details$cond_type == "aa_change" && !is.na(cond_details$AA_ref)) {
          if(a_details$start_AA_ref == cond_details$AA_ref) {
            # Current insertion variant start position matches a defined variant, but we can't be sure we have an exact match to the defined AA change
            manual_review <- TRUE
            next
          }
        }
      } else if(a_details$mut_type == "del" & any((a_details$start_pos:a_details$end_pos) %in% (cond_details$start_pos:cond_details$end_pos))) {
        # Current variant is a deletion and deleted positions overlap defined range
        if(cond_details$cond_type == "aa_range") {
          # Current variant in defined range
          whitelist <- TRUE
          break
        } else if(cond_details$cond_type == "aa_change") {
          # Current deletion variant affects this position, but we can't be sure we have an exact match to the defined AA change
          manual_review <- TRUE
          next
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

match_ns <- function(x) {
  whitelist <- FALSE
  manual_review <- FALSE
  putative_whitelist <- FALSE
  putative_manual_review <- FALSE
  if(grepl("(nonsynonymous|nonframeshift)", x[[1]], perl = TRUE) && x[[4]] != "nan" && x[[5]] != "nan") {
    conditions <- strsplit(x[[2]], ";")[[1]]
    conditions_put <- strsplit(x[[3]], ";")[[1]]
    aa_changes <- strsplit(x[[4]], ",")[[1]]
    aa_exons <- as.integer(strsplit(x[[5]], ",")[[1]])
    match_conditions <- match_ns_conditions(conditions, aa_changes, aa_exons)
    match_put_conditions <- match_ns_conditions(conditions_put, aa_changes, aa_exons)
    whitelist <- match_conditions$whitelist
    manual_review <- match_conditions$manual_review
    putative_whitelist <- match_put_conditions$whitelist
    putative_manual_review <- match_put_conditions$manual_review
  }
  return(c(whitelist, manual_review, putative_whitelist, putative_manual_review))
}

match_fs_conditions <- function(conditions, aa_changes, aa_exons) {
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
    }
    cond_details <- parse_cond(cond)
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
      a_details <- parse_aa_change(a)
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

match_fs <- function(x) {
  whitelist <- FALSE
  manual_review <- FALSE
  putative_whitelist <- FALSE
  putative_manual_review <- FALSE
  if(grepl("^frameshift", x[[1]], perl = TRUE) && x[[4]] != "nan" && x[[5]] != "nan") {
    conditions <- strsplit(x[[2]], ";")[[1]]
    conditions_put <- strsplit(x[[3]], ";")[[1]]
    aa_changes <- strsplit(x[[4]], ",")[[1]]
    aa_exons <- as.integer(strsplit(x[[5]], ",")[[1]])
    match_conditions <- match_fs_conditions(conditions, aa_changes, aa_exons)
    match_put_conditions <- match_fs_conditions(conditions_put, aa_changes, aa_exons)
    whitelist <- match_conditions$whitelist
    manual_review <- match_conditions$manual_review
    putative_whitelist <- match_put_conditions$whitelist
    putative_manual_review <- match_put_conditions$manual_review
  }
  return(c(whitelist, manual_review, putative_whitelist, putative_manual_review))
}

match_sg_conditions <- function(conditions, aa_changes, aa_exons) {
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
    }
    cond_details <- parse_cond(cond)
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
      a_details <- parse_aa_change(a)
      if(all(is.na(a_details))) {
        next
      }
      if(a_details$mut_type == "aa_change" & (a_details$start_pos %in% cond_details$start_pos:cond_details$end_pos)) {
        # Current variant is a single AA substitution
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
      } else if(a_details$mut_type == "ins" & (a_details$start_pos %in% cond_details$start_pos:cond_details$end_pos)) {
        # Current variant is an insertion
        if(cond_details$cond_type == "aa_range") {
          # Current variant in defined range
          whitelist <- TRUE
          break
        } else if(cond_details$cond_type == "aa_change" && !is.na(cond_details$AA_ref)) {
          if(a_details$start_AA_ref == cond_details$AA_ref) {
            # Current insertion variant start position matches a defined variant, but we can't be sure we have an exact match to the defined AA change
            manual_review <- TRUE
            next
          }
        }
      } else if(a_details$mut_type == "del" & any((a_details$start_pos:a_details$end_pos) %in% (cond_details$start_pos:cond_details$end_pos))) {
        # Current variant is a deletion and deleted positions overlap defined range
        if(cond_details$cond_type == "aa_range") {
          # Current variant in defined range
          whitelist <- TRUE
          break
        } else if(cond_details$cond_type == "aa_change") {
          # Current deletion variant affects this position, but we can't be sure we have an exact match to the defined AA change
          manual_review <- TRUE
          next
        }
      } else if(a_details$mut_type == "fs" & (a_details$start_pos %in% cond_details$start_pos:cond_details$end_pos)) {
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

match_sg <- function(x) {
  whitelist <- FALSE
  manual_review <- FALSE
  putative_whitelist <- FALSE
  putative_manual_review <- FALSE
  if(grepl("stopgain", x[[1]], perl = TRUE) && x[[4]] != "nan" && x[[5]] != "nan") {
    conditions <- strsplit(x[[2]], ";")[[1]]
    conditions_put <- strsplit(x[[3]], ";")[[1]]
    aa_changes <- strsplit(x[[4]], ",")[[1]]
    aa_exons <- as.integer(strsplit(x[[5]], ",")[[1]])
    match_conditions <- match_sg_conditions(conditions, aa_changes, aa_exons)
    match_put_conditions <- match_sg_conditions(conditions_put, aa_changes, aa_exons)
    whitelist <- match_conditions$whitelist
    manual_review <- match_conditions$manual_review
    putative_whitelist <- match_put_conditions$whitelist
    putative_manual_review <- match_put_conditions$manual_review
  }
  return(c(whitelist, manual_review, putative_whitelist, putative_manual_review))
}

match_sp_conditions <- function(conditions, transcript_changes, transcript_exons) {
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
      }
      cond_details <- parse_cond(cond)
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

match_sp <- function(x) {
  whitelist <- FALSE
  manual_review <- FALSE
  putative_whitelist <- FALSE
  putative_manual_review <- FALSE
  if(grepl("splicing", x[[1]], perl = TRUE) && x[[4]] != "nan" && x[[5]] != "nan") {
    conditions <- strsplit(x[[2]], ";")[[1]]
    conditions_put <- strsplit(x[[3]], ";")[[1]]
    transcript_changes <- strsplit(x[[4]], ";")[[1]]
    transcript_exons <- as.integer(strsplit(x[[5]], ";")[[1]])
    match_conditions <- match_sp_conditions(conditions, transcript_changes, transcript_exons)
    match_put_conditions <- match_sp_conditions(conditions_put, transcript_changes, transcript_exons)
    whitelist <- match_conditions$whitelist
    manual_review <- match_conditions$manual_review
    putative_whitelist <- match_put_conditions$whitelist
    putative_manual_review <- match_put_conditions$manual_review
  }
  return(c(whitelist, manual_review, putative_whitelist, putative_manual_review))
}

# Matching non-synonymous mutations (SNV, non-frameshift INDELs)
wl_ns <- apply(vars_g_chip_func_filtered[, c("ExonicFunc.refGene", "Nonsynonymous", "Putative_Nonsynonymous", "AAChange.protChange", "AAChange.exon")], 1, match_ns)

# Matching frameshift mutations
wl_fs <- apply(vars_g_chip_func_filtered[, c("ExonicFunc.refGene", "Frameshift", "Putative_Frameshift", "AAChange.protChange", "AAChange.exon")], 1, match_fs)

# Matching stop-gain mutations
wl_sg <- apply(vars_g_chip_func_filtered[, c("ExonicFunc.refGene", "Stop_Gain", "Putative_Stop_Gain", "AAChange.protChange", "AAChange.exon")], 1, match_sg)

# Matching splicing mutations
wl_sp <- apply(vars_g_chip_func_filtered[, c("Func.refGene", "Splicing", "Putative_Splicing", "GeneDetail.transcript", "GeneDetail.exon")], 1, match_sp)

# Add whitelist/manual review columns to data frame
wl_ns_df <- data.frame(t(wl_ns))
colnames(wl_ns_df) <- paste(c("Whitelist", "Manual_Review", "Putative_Whitelist", "Putative_Manual_Review"), "Nonsynonymous", sep = "_")
wl_fs_df <- data.frame(t(wl_fs))
colnames(wl_fs_df) <- paste(c("Whitelist", "Manual_Review", "Putative_Whitelist", "Putative_Manual_Review"), "Frameshift", sep = "_")
wl_sg_df <- data.frame(t(wl_sg))
colnames(wl_sg_df) <- paste(c("Whitelist", "Manual_Review", "Putative_Whitelist", "Putative_Manual_Review"), "Stop_Gain", sep = "_")
wl_sp_df <- data.frame(t(wl_sp))
colnames(wl_sp_df) <- paste(c("Whitelist", "Manual_Review", "Putative_Whitelist", "Putative_Manual_Review"), "Splicing", sep = "_")

# Overall whitelist/manual review
wl <- wl_ns_df[, 1] | wl_fs_df[, 1] | wl_sg_df[, 1] | wl_sp_df[, 1]
mr <- (!wl) & (wl_ns_df[, 2] | wl_fs_df[, 2] | wl_sg_df[, 2] | wl_sp_df[, 2])
pwl <- wl_ns_df[, 3] | wl_fs_df[, 3] | wl_sg_df[, 3] | wl_sp_df[, 3]
pmr <- (!pwl) & (wl_ns_df[, 4] | wl_fs_df[, 4] | wl_sg_df[, 4] | wl_sp_df[, 4])
overall_wl <- data.frame(Whitelist = wl, Manual_Review = mr, Putative_Whitelist = pwl, Putative_Manual_Review = pmr)

vars_g_chip_func_filtered_wl <- cbind(vars_g_chip_func_filtered, wl_ns_df, wl_fs_df, wl_sg_df, wl_sp_df, overall_wl)

write.csv(vars_g_chip_func_filtered_wl, paste(sample_id, ".all_variants.csv", sep = ""), row.names = FALSE)
write.csv(vars_g_chip_func_filtered_wl[vars_g_chip_func_filtered_wl$Whitelist,], paste(sample_id, ".chip_wl_variants.csv", sep = ""), row.names = FALSE)
write.csv(vars_g_chip_func_filtered_wl[vars_g_chip_func_filtered_wl$Manual_Review,], paste(sample_id, ".chip_manual_review_variants.csv", sep = ""), row.names = FALSE)
write.csv(vars_g_chip_func_filtered_wl[vars_g_chip_func_filtered_wl$Putative_Whitelist,], paste(sample_id, ".chip_putative_wl_variants.csv", sep = ""), row.names = FALSE)
write.csv(vars_g_chip_func_filtered_wl[vars_g_chip_func_filtered_wl$Putative_Manual_Review,], paste(sample_id, ".chip_putative_manual_review_variants.csv", sep = ""), row.names = FALSE)
