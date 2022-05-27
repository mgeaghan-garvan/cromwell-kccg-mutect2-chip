# library(data.table, quietly = T)
library(tidyr)
library(stringr)

# Source R scripts
source("./import/gnomad.R")
source("./import/hard_filter.R")
source("./import/homopolymers.R")
source("./import/parse_annovar.R")
source("./import/match_mutation.R")

source("./import/update_whitelist.R")
source("./import/protein_terminal.R")

# source("./import/extract_var_details.R")
# source("./import/parse_chip_defs.R")
# source("./import/match_nonsynonymous.R")
# source("./import/match_frameshift.R")
# source("./import/match_stopgain.R")
# source("./import/match_splicing.R")
# source("./import/apply_chip_filters.R")
# source("./import/apply_putative_filter.R")


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


# ================== #
# Apply hard filters #
# ================== #

# Filter by AD, DP, AF, F1R2/F2R1 VCF fields and gnomAD frequency
vars_hf <- apply_hard_filters(vars, tumor_sample_name)


# ================================================ #
# Apply filters to variants in homopolymer regions #
# ================================================ #

# NOTE: For now, we'll just be doing INDELs (HP_SNVS = FALSE), but in the future we should make this configurable via CLI arguments
vars_hf_hp <- apply_homopolymer_filter(vars_hf, fasta_file, HP_SNVS = FALSE)


# ==================== #
# PARSE ANNOVAR OUTPUT #
# ==================== #

# Parse annovar output, split rows into one per transcript, and extract transcript accession IDs
vars_ann <- parse_annovar(vars_hf_hp)

# Filter for CHIP transcripts
refseq_ensembl_suffix = "ensGene"  # change this to "refGene" if you want to use RefSeq IDs
refseq_ensembl_chip_accession_column = "ensembl_accession"  # change this to "refseq_accession" if you want to use RefSeq IDs
aachange <- paste("AAChange", refseq_ensembl_suffix, sep = ".")
genedetail <- paste("GeneDetail", refseq_ensembl_suffix, sep = ".")
transcript <- paste("Transcript", refseq_ensembl_suffix, sep = ".")

vars_chip <- vars_ann[vars_ann[[transcript]] %in% chip_vars[[refseq_ensembl_chip_accession_column]],]

# Merge filtered variants and CHIP gene annotation data frames
vars_chip <- merge(vars_chip, chip_vars, by.x = transcript, by.y = refseq_ensembl_chip_accession_column)


# ============================================================================== #
# Filter rows if chip mutation definitions match AAChange/GeneDetail column info #
# ============================================================================== #

vars_chip_filtered <- match_mut_def(vars_chip)


# ==================================================== #
# Apply stricter hard filter to putative CHIP variants #
# ==================================================== #


# ========== OLD CODE ========== #


vars_chip$AAChange.protChange <- unlist(lapply(vars_chip[, c("AAChange.transcript")], extractNonsyn))

vars_chip$AAChange.exon <- unlist(lapply(vars_chip[, c("AAChange.transcript")], extractNonsynExon))

vars_chip$GeneDetail.exon <- unlist(lapply(vars_chip[, c("GeneDetail.transcript")], extractSpliceExons ))

# Remove variants with no matching affected transcript
vars_chip_filtered <- vars_chip[!(
  vars_chip$GeneDetail.transcript == "nan" &
    vars_chip$AAChange.transcript == "nan" &
    vars_chip$AAChange.protChange == "nan" &
    vars_chip$AAChange.exon == "nan" &
    vars_chip$GeneDetail.exon == "nan" &
    !(vars_chip$Nonsynonymous %in% c("any", "c_term")) &
    !(vars_chip$Frameshift %in% c("any", "c_term")) &
    !(vars_chip$Stop_Gain %in% c("any", "c_term")) &
    !(vars_chip$Splicing %in% c("any", "c_term")) &
    !(vars_chip$Putative_Nonsynonymous %in% c("any", "c_term")) &
    !(vars_chip$Putative_Frameshift %in% c("any", "c_term")) &
    !(vars_chip$Putative_Stop_Gain %in% c("any", "c_term")) &
    !(vars_chip$Putative_Splicing %in% c("any", "c_term"))
), ]


# ============================================================== #
# Get position along protein - add N- and C-terminal 10% columns #
# ============================================================== #
vars_chip_filtered <- get_protein_terminal_status(vars_chip_filtered, transcript_prot_file, parse_aa_change)


# ==================================================================== #
# Apply CHIP variant definition filters and generate initial whitelist #
# ==================================================================== #

whitelist <- apply_chip_filters(vars_chip_filtered, match_ns, match_fs, match_sg, match_sp, parse_cond, parse_aa_change)


# ==================================================== #
# Apply stricter hard filter to putative CHIP variants #
# ==================================================== #

vars_chip_filtered <- apply_putative_filter(vars_chip_filtered, whitelist)

# Update the overall whitelist filter columns to include the putative variant filter information
whitelist$overall_wl <- update_whitelist(whitelist$overall_wl, vars_chip_filtered$COMBINED_FILTER)


# ================== #
# Finalise dataframe #
# ================== #

# Merge whitelist columns into main dataframe
vars_chip_filtered_wl <- cbind(vars_chip_filtered, whitelist$wl_ns_df, whitelist$wl_fs_df, whitelist$wl_sg_df, whitelist$wl_sp_df, whitelist$overall_wl)

# Perform final removal of duplicate entries
vars_chip_filtered_wl <- unique(vars_chip_filtered_wl)


# ============= #
# Write to file #
# ============= #

write.csv(vars_chip_filtered_wl, paste(sample_id, ".all_variants.csv", sep = ""), row.names = FALSE)
write.csv(vars_chip_filtered_wl[vars_chip_filtered_wl$Whitelist,], paste(sample_id, ".chip_wl_variants.csv", sep = ""), row.names = FALSE)
write.csv(vars_chip_filtered_wl[vars_chip_filtered_wl$Manual_Review,], paste(sample_id, ".chip_manual_review_variants.csv", sep = ""), row.names = FALSE)
write.csv(vars_chip_filtered_wl[vars_chip_filtered_wl$Putative_Whitelist,], paste(sample_id, ".chip_putative_wl_variants.csv", sep = ""), row.names = FALSE)
write.csv(vars_chip_filtered_wl[vars_chip_filtered_wl$Putative_Manual_Review,], paste(sample_id, ".chip_putative_manual_review_variants.csv", sep = ""), row.names = FALSE)
