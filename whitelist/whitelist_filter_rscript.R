# library(data.table, quietly = T)
library(tidyr)
library(stringr)

# Source R scripts
source("./import/gnomad.R")
source("./import/hard_filter.R")
source("./import/homopolymers.R")
source("./import/extract_var_details.R")
source("./import/parse_chip_defs.R")
source("./import/match_nonsynonymous.R")
source("./import/match_frameshift.R")
source("./import/match_stopgain.R")
source("./import/match_splicing.R")
source("./import/apply_chip_filters.R")
source("./import/apply_putative_filter.R")
source("./import/update_whitelist.R")
source("./import/protein_terminal.R")


# ========================== #
# Get command-line arguments #
# ========================== #

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9) {
  stop(paste(
    "Incorrect number of arguments!",
    "Usage: Rscript whitelist_filter_rscript.R <ANNOVAR OUTPUT TABLE> <ANNOVAR OUTPUT VCF> <TUMOR SAMPLE NAME> <GNOMAD SOURCE> <GNOMAD SUBPOPULATION CODE> <TREAT MISSING AF AS RARE> <CHIP DEFINITION FILE> <TRANSCRIPT PROTEIN LENGTHS FILE> <FASTA REFERENCE FILE>",
    "    ANNOVAR OUTPUT TABLE/VCF:         output txt and vcf files from Annovar.",
    "    TUMOR SAMPLE NAME:                name of the tumor sample as recorded in the VCF file's column header line.",
    "    GNOMAD SOURCE:                    one of 'exome', 'genome', 'exome,genome', or 'genome,exome' (must be in the same order as when running annovar)",
    "    GNOMAD SUBPOPULATION CODE:        one of 'AF' (all), 'AF_afr', 'AF_sas', 'AF_amr', 'AF_eas', 'AF_nfe', 'AF_fin', 'AF_asj'",
    "    TREAT MISSING AF AS RARE:         'TRUE' = variants not annotated in gnomAD are assumed to be rare and are given an allele frequency of 0; 'FALSE' = variants not annotated in gnomAD will not pass the gnomAD hard filter.",
    "    CHIP DEFINITION FILE:             csv file containing CHIP variant definitions",
    "    TRANSCRIPT PROTEIN LENGTHS FILE:  a file containing three columns: RefSeq transcript ID, HGNC gene symbol, protein length",
    "    FASTA REFERENCE FILE:             the FASTA reference to which the samples have been aligned.",
    sep = "\n"))
}
annovar_text_out <- args[1]
annovar_vcf_out <- args[2]
tumor_sample_name <- args[3]
gnomad_source <- args[4]
gnomad_pop <- args[5]
treat_missing_as_rare <- args[6]
chip_def_file <- args[7]
transcript_prot_file <- args[8]
fasta_file <- args[9]


# ============================ #
# Check command-line arguments #
# ============================ #

if (!file.exists(annovar_text_out)) {
  stop("Input annovar table file does not exist.")
}
if (!file.exists(annovar_vcf_out)) {
  stop("Input annovar vcf file does not exist.")
}
if (!file.exists(chip_def_file)) {
  stop("Chip variant definition file does not exist.")
}
if (!file.exists(transcript_prot_file)) {
  stop("Protein length file does not exist.")
}
if (!file.exists(fasta_file)) {
  stop("FASTA reference file does not exist.")
}
if (!(gnomad_source %in% c("exome", "genome", "exome,genome", "genome,exome"))) {
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


# ========= #
# Load data #
# ========= #

# Load CHIP variant/gene lists
chip_vars <- read.csv(chip_def_file)

# Define sample ID
sample_id <- gsub(annovar_text_out_regex, "", annovar_text_out, perl = TRUE)
sample_id <- gsub("^.*\\/([^\\/]+)$", "\\1", sample_id, perl = TRUE)

# Load annovar variant annotations
vars <- read.delim(annovar_text_out, sep = "\t", header = TRUE)
vars$Sample <- sample_id
vars$TumorSample <- tumor_sample_name

# Load annovar vcf file
vcf <- scan(annovar_vcf_out, character(), sep = "\n")
vcf_header <- grep("^#CHROM", vcf, perl = TRUE, value = TRUE)
vcf_header <- gsub("^#", "", vcf_header, perl = TRUE)
vcf_header <- strsplit(vcf_header, "\t")[[1]]


# ============== #
# Rename columns #
# ============== #

# Rename Otherinfo columns
# NOTE: These are hard-coded to match the output of annovar as of 2021-11-22
#       When converting from vcf to avinput format, the following command is used by table_annovar.pl:
#       convert2annovar.pl -includeinfo -allsample -withfreq -format vcf4 INPUT.vcf > OUTPUT.avinput
#       This produces 3 "Otherinfo" columns corresponding to total allele frequency across all samples, quality score, and read deapth
#       The following 9+N columns (for N samples) correspond to the input vcf file's columns
otherinfo_cols <- c("ANNOVAR_alt_af", "ANNOVAR_qual", "ANNOVAR_alt_ad", vcf_header)
colnames(vars)[grepl("Otherinfo", colnames(vars), fixed = TRUE)] <- otherinfo_cols

# Rename gnomad columns and get gnomad AF
vars <- rename_gnomad_col(vars, gnomad_source)
vars <- get_gnomad_af(vars, gnomad_source, gnomad_pop, treat_missing_as_rare)

# ============================================= #
# Define regex strings for mutation definitions #
# ============================================= #

regex_list <- list(
  all = "^all$",
  c_terminal = "^c_term$",
  exon_any = "^exon\\d+$",
  exon_aa_insertion_specific = "^exon\\d+\\[ins[A-Z]\\]$",
  range = "^\\d+\\-\\d+$",
  substitution_specific = "^[A-Z]\\d+[A-Z]$",
  substitution_any = "^[A-Z]\\d+$",
  deletion_single = "^[A-Z]\\d+del$",
  deletion_range = "^[A-Z]\\d+_[A-Z]\\d+del$",
  insertion_specific = "^[A-Z]\\d+_[A-Z]\\d+ins[A-Z]+$",
  insertion_any = "^[A-Z]\\d+_[A-Z]\\d+ins$",
  deletion_insertion = "^[A-Z]\\d+_[A-Z]\\d+delins[A-Z]+$",
  coding_range = "^c\\.\\d+-\\d+$",
  coding_insertion_any = "^c\\.\\d+_\\d+ins$",
  frameshift_any = "^[A-Z]\\d+fs$",
  stop_gain = "^[A-Z]\\d+\\*$"
)


# ================== #
# Apply hard filters #
# ================== #

# Filter by AD, DP, AF, F1R2/F2R1 VCF fields and gnomAD frequency
vars <- apply_hard_filters(vars, tumor_sample_name)


# ================================================ #
# Apply filters to variants in homopolymer regions #
# ================================================ #

# NOTE: For now, we'll just be doing INDELs (HP_SNVS = FALSE), but in the future we should make this configurable via CLI arguments
vars <- apply_homopolymer_filter(vars, fasta_file, HP_SNVS = FALSE)


# ==================== #
# PARSE ANNOVAR OUTPUT #
# ==================== #

# Split rows by RefSeq IDs in both GeneDetail.refGene (noncoding variant annotations, e.g. splicing variants) and AAChange.refGene (coding variant annotations)
vars_g <- as.data.frame(separate_rows(vars, GeneDetail.refGene, sep = ";"))
vars_g <- as.data.frame(separate_rows(vars_g, AAChange.refGene, sep = ","))

# Filter variants for exonic or splicing variants
vars_g_func <- vars_g[
  (
    grepl("exonic", vars_g$Func.refGene, fixed = T) |
      grepl("splicing", vars_g$Func.refGene, fixed = T)
  ),
]

# Drop potential duplicates
vars_g_func <- unique(vars_g_func)

# Construct a column with the RefSeq IDs from both GeneDetail.refGene and AAChange.refGene, and split rows again on this column
vars_g_func_acc <- unlist(apply(vars_g_func[c("AAChange.refGene", "GeneDetail.refGene")], 1, function(x) {
  aa <- as.character(x[[1]])
  gd <- as.character(x[[2]])
  aagd <- paste(aa, gd, sep = ":")
  ann <- strsplit(aagd, ":")[[1]]
  acc <- grep("^N[MR]_\\d+(\\.\\d+)?$", ann, perl = TRUE, value = TRUE)
  acc <- gsub("^(N[MR]_\\d+)(\\.\\d+)?$", "\\1", acc, perl = TRUE)
  return(paste(acc, collapse = ";"))
}))
vars_g_func$Transcript.refGene <- vars_g_func_acc
vars_g_func <- as.data.frame(separate_rows(vars_g_func, Transcript.refGene, sep = ";"))
acc_func <- function(x) {
  acc <- as.character(x[[1]])
  ann <- as.character(x[[2]])
  ret <- grep(paste("(^|:)", acc, "(:|$)", sep = ""), ann, perl = TRUE, value = TRUE)
  if (length(ret) != 1) {
    return("")
  } else {
    return(ret)
  }
}
vars_g_func$AAChange.refGene <- unlist(apply(vars_g_func[c("Transcript.refGene", "AAChange.refGene")], 1, acc_func))
vars_g_func$GeneDetail.refGene <- unlist(apply(vars_g_func[c("Transcript.refGene", "GeneDetail.refGene")], 1, acc_func))

# Drop potential duplicates and rows without a RefSeq ID
vars_g_func <- unique(vars_g_func[vars_g_func$Transcript.refGene != "",])


# Filter for CHIP transcripts
vars_g_chip_func <- vars_g_func[vars_g_func$Transcript.refGene %in% chip_vars$refseq_accession,]


# Merge filtered variants and CHIP gene annotation data frames
vars_g_chip_func <- merge(vars_g_chip_func, chip_vars, by.x = "Transcript.refGene", by.y = "refseq_accession")


# =============================================== #
# Extract variant information from Annovar output #
# =============================================== #

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


# ============================================================== #
# Get position along protein - add N- and C-terminal 10% columns #
# ============================================================== #
vars_g_chip_func_filtered <- get_protein_terminal_status(vars_g_chip_func_filtered, transcript_prot_file, parse_aa_change)


# ==================================================================== #
# Apply CHIP variant definition filters and generate initial whitelist #
# ==================================================================== #

whitelist <- apply_chip_filters(vars_g_chip_func_filtered, match_ns, match_fs, match_sg, match_sp, parse_cond, parse_aa_change)


# ==================================================== #
# Apply stricter hard filter to putative CHIP variants #
# ==================================================== #

vars_g_chip_func_filtered <- apply_putative_filter(vars_g_chip_func_filtered, whitelist)

# Update the overall whitelist filter columns to include the putative variant filter information
whitelist$overall_wl <- update_whitelist(whitelist$overall_wl, vars_g_chip_func_filtered$COMBINED_FILTER)


# ================== #
# Finalise dataframe #
# ================== #

# Merge whitelist columns into main dataframe
vars_g_chip_func_filtered_wl <- cbind(vars_g_chip_func_filtered, whitelist$wl_ns_df, whitelist$wl_fs_df, whitelist$wl_sg_df, whitelist$wl_sp_df, whitelist$overall_wl)

# Perform final removal of duplicate entries
vars_g_chip_func_filtered_wl <- unique(vars_g_chip_func_filtered_wl)


# ============= #
# Write to file #
# ============= #

write.csv(vars_g_chip_func_filtered_wl, paste(sample_id, ".all_variants.csv", sep = ""), row.names = FALSE)
write.csv(vars_g_chip_func_filtered_wl[vars_g_chip_func_filtered_wl$Whitelist,], paste(sample_id, ".chip_wl_variants.csv", sep = ""), row.names = FALSE)
write.csv(vars_g_chip_func_filtered_wl[vars_g_chip_func_filtered_wl$Manual_Review,], paste(sample_id, ".chip_manual_review_variants.csv", sep = ""), row.names = FALSE)
write.csv(vars_g_chip_func_filtered_wl[vars_g_chip_func_filtered_wl$Putative_Whitelist,], paste(sample_id, ".chip_putative_wl_variants.csv", sep = ""), row.names = FALSE)
write.csv(vars_g_chip_func_filtered_wl[vars_g_chip_func_filtered_wl$Putative_Manual_Review,], paste(sample_id, ".chip_putative_manual_review_variants.csv", sep = ""), row.names = FALSE)
