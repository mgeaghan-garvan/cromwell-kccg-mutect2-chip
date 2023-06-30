# library(data.table, quietly = T)
library(tidyr)
library(stringr)
library(dplyr)

# Source R scripts
source("./import/gnomad.R")
source("./import/hard_filter.R")
source("./import/homopolymers.R")
source("./import/somaticism_filter.R")
source("./import/parse_annovar.R")
source("./import/match_mutation.R")
source("./import/apply_putative_filter.R")


# ========================== #
# Get command-line arguments #
# ========================== #

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 12) {
  stop(paste(
    "Incorrect number of arguments!",
    "Usage: Rscript whitelist_filter_rscript.R <INPUT_VCF> <ANNOVAR OUTPUT TABLE> <ANNOVAR VARIANT FUNCTION TABLE> <ANNOVAR VARIANT EXONIC FUNCTION TABLE> <ENSEMBL REFSEQ> <TUMOR SAMPLE NAME> <GNOMAD SOURCE> <GNOMAD SUBPOPULATION CODE> <TREAT MISSING AF AS RARE> <CHIP DEFINITION FILE> <FASTA REFERENCE FILE> <SOMATICISM FILTER TRANSCRIPTS>",
    "    INPUT_VCF:                               VCF file to annotate for CHIP.",
    "    ANNOVAR OUTPUT TABLE:                    output txt files from Annovar.",
    "    ANNOVAR VARIANT FUNCTION TABLE:          variant function output file from Annovar.",
    "    ANNOVAR VARIANT EXONIC FUNCTION TABLE:   variant exonic function output file from Annovar.",
    "    ENSEMBL REFSEQ:                          either 'ensembl' or 'refseq', specifying which annotation to use when matching against CHIP definitions.",
    "    TUMOR SAMPLE NAME:                       name of the tumor sample as recorded in the VCF file's column header line.",
    "    GNOMAD SOURCE:                           one of 'exome', 'genome', 'exome,genome', or 'genome,exome' (must be in the same order as when running annovar)",
    "    GNOMAD SUBPOPULATION CODE:               one of 'AF' (all), 'AF_afr', 'AF_sas', 'AF_amr', 'AF_eas', 'AF_nfe', 'AF_fin', 'AF_asj'",
    "    TREAT MISSING AF AS RARE:                'TRUE' = variants not annotated in gnomAD are assumed to be rare and are given an allele frequency of 0; 'FALSE' = variants not annotated in gnomAD will not pass the gnomAD hard filter.",
    "    CHIP DEFINITION FILE:                    csv file containing CHIP variant definitions",
    "    FASTA REFERENCE FILE:                    the FASTA reference to which the samples have been aligned.",
    "    SOMATICISM FILTER TRANSCRIPTS:           text file containing transcript IDs (one per line) to be subjected to the somaticism filter (to remove likely germline mutations).",
    sep = "\n"))
}
input_vcf <- args[1]
annovar_text_out <- args[2]
annovar_variant_func_out <- args[3]
annovar_variant_exonic_func_out <- args[4]
ensembl_refseq <- args[5]
tumor_sample_name <- args[6]
gnomad_source <- args[7]
gnomad_pop <- args[8]
treat_missing_as_rare <- args[9]
chip_def_file <- args[10]
fasta_file <- args[11]
somaticism_file <- args[12]


# ============================ #
# Check command-line arguments #
# ============================ #

if (!file.exists(input_vcf)) {
  stop("Input VCF file does not exist.")
}
if (!file.exists(annovar_text_out)) {
  stop("Input annovar table file does not exist.")
}
if (!file.exists(annovar_variant_func_out)) {
  stop("Input annovar variant function file does not exist.")
}
if (!file.exists(annovar_variant_exonic_func_out)) {
  stop("Input annovar variant exonic function file does not exist.")
}
if (!(ensembl_refseq %in% c("ensembl", "refseq"))) {
  stop("Must specify either 'ensembl' or 'refseq' annotation to use.")
}
if (!file.exists(chip_def_file)) {
  stop("Chip variant definition file does not exist.")
}
if (!file.exists(fasta_file)) {
  stop("FASTA reference file does not exist.")
}
if (!file.exists(somaticism_file)) {
  stop("Somaticism filter transcripts file does not exist.")
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


# ========= #
# Load data #
# ========= #

# Load CHIP variant/gene lists
chip_vars <- read.csv(chip_def_file)

# Remove version numbers from RefSeq and Ensembl accession IDs if present
# i.e. NM12345.1 -> NM12345
chip_vars[["refseq_accession"]] <- gsub("\\.\\d+$", "", chip_vars[["refseq_accession"]])
chip_vars[["ensembl_accession"]] <- gsub("\\.\\d+$", "", chip_vars[["ensembl_accession"]])

# Define sample ID
sample_id <- gsub(annovar_text_out_regex, "", annovar_text_out, perl = TRUE)
sample_id <- gsub("^.*\\/([^\\/]+)$", "\\1", sample_id, perl = TRUE)

# Load annovar variant annotations
vars <- read.delim(annovar_text_out, sep = "\t", header = TRUE)
vars$Sample <- sample_id
vars$TumorSample <- tumor_sample_name

# Load input vcf file
vcf_header = ""
con <- file(input_vcf, "r")
while (TRUE) {
  line <- readLines(con, n = 1)
  if (length(line) == 0) {
    break
  } else if (grepl("^#CHROM", line, perl = TRUE)) {
    vcf_header <- strsplit(gsub("^#", "", line, perl = TRUE), "\t")[[1]]
    break
  }
}

# vcf <- scan(input_vcf, character(), sep = "\n")
# vcf_header <- grep("^#CHROM", vcf, perl = TRUE, value = TRUE)
# vcf_header <- gsub("^#", "", vcf_header, perl = TRUE)
# vcf_header <- strsplit(vcf_header, "\t")[[1]]
# vcf_whole_header <- grep("^#", vcf, perl = TRUE, value = TRUE)
# vcf_body <- grep("^[^#]", vcf, perl = TRUE, value = TRUE)
# vcf_body_df <- as.data.frame(do.call(rbind, strsplit(vcf_body, split = "\t", fixed = TRUE)))
# colnames(vcf_body_df) <- vcf_header


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

# Index variants for later annotation of VCF
vars_annovar_id_str <- paste(vars$Chr, vars$Start, vars$End, vars$Ref, vars$Alt, sep = ":")
vars_vcf_id_str <- paste(vars$CHR, vars$POS, vars$REF, vars$ALT, sep = ":")
vars_annovar_vcf_ids <- data.frame(annovar_id = vars_annovar_id_str, vcf_id = vars_vcf_id_str, idx = 1, row.names = vars_annovar_id_str)
tmp_seen <- list()
for (idx in 1:length(vars_annovar_id_str)) {
  annovar_id <- vars_annovar_id_str[idx]
  vcf_id <- vars_vcf_id_str[idx]
  if (is.null(tmp_seen[[vcf_id]])) {
    tmp_seen[[vcf_id]] <- TRUE
  } else {
    vars_annovar_vcf_ids[annovar_id, 'idx'] <- vars_annovar_vcf_ids[annovar_id, 'idx'] + 1
  }
}
rm(tmp_seen)


# Rename gnomad columns and get gnomad AF
vars <- rename_gnomad_col(vars, gnomad_source)
vars <- get_gnomad_af(vars, gnomad_source, gnomad_pop, treat_missing_as_rare)


# ================== #
# Apply hard filters #
# ================== #

# Filter by AD, DP, AF, F1R2/F2R1 VCF fields and gnomAD frequency
vars_hf <- apply_hard_filters(vars, tumor_sample_name)


# =========================================== #
# Select between using RefSeq and Ensembl IDs #
# =========================================== #

# ensembl_refseq <- "refseq"  # Use Ensembl IDs by default. Change this to "refseq" if you want to use RefSeq IDs
ensGene <- TRUE
refGene <- FALSE
refseq_ensembl_suffix <- "ensGene"
refseq_ensembl_chip_accession_column <- "ensembl_accession"
refseq_ensembl_other_chip_accession_column <- "refseq_accession"
refseq_ensembl_variant_func_regex_sub <- "(ENST\\d+)([^\\d]|$)"
refseq_ensembl_variant_func_regex_match <- "^ENST\\d+$"
if (ensembl_refseq == "refseq") {
  ensGene <- FALSE
  refGene <- TRUE
  refseq_ensembl_suffix <- "refGene"
  refseq_ensembl_chip_accession_column <- "refseq_accession"
  refseq_ensembl_other_chip_accession_column <- "ensembl_accession"
  refseq_ensembl_variant_func_regex_sub <- "(NM_\\d+)([^\\d]|$)"
  refseq_ensembl_variant_func_regex_match <- "^NM_\\d+$"
}

aachange <- paste("AAChange", refseq_ensembl_suffix, sep = ".")
genedetail <- paste("GeneDetail", refseq_ensembl_suffix, sep = ".")
transcript <- paste("Transcript", refseq_ensembl_suffix, sep = ".")
exonic_func <- paste("ExonicFunc", refseq_ensembl_suffix, sep = ".")
func <- paste("Func", refseq_ensembl_suffix, sep = ".")


# ================================================ #
# Apply filters to variants in homopolymer regions #
# ================================================ #

vars_hf_hp <- apply_homopolymer_indel_filter(vars_hf, fasta_file)

# ==================== #
# PARSE ANNOVAR OUTPUT #
# ==================== #

# Load annovar variant function annotations
vars_variant_func <- read.delim(annovar_variant_func_out, sep = "\t", header = FALSE)
vars_variant_func <- unique(vars_variant_func[1:7])
colnames(vars_variant_func) <- c("Func", "Transcript", "Chr", "Start", "End", "Ref", "Alt")
vars_variant_func$Transcript <- gsub(refseq_ensembl_variant_func_regex_sub, ",\\1,", vars_variant_func$Transcript, perl = TRUE)
vars_variant_func$Transcript <- gsub("^,", "", vars_variant_func$Transcript, perl = TRUE)
vars_variant_func$Transcript <- gsub(",,+", ",", vars_variant_func$Transcript, perl = TRUE)
vars_variant_func <- (vars_variant_func %>% separate_rows(Transcript, sep = ","))
vars_variant_func <- vars_variant_func[grepl(refseq_ensembl_variant_func_regex_match, vars_variant_func$Transcript),]
vars_variant_func <- unique(vars_variant_func)

# Load annovar variant exonic function annotations
vars_variant_exonic_func <- read.delim(annovar_variant_exonic_func_out, sep = "\t", header = FALSE)
vars_variant_exonic_func <- unique(vars_variant_exonic_func[2:8])
colnames(vars_variant_exonic_func) <- c("ExonicFunc", "AAChange", "Chr", "Start", "End", "Ref", "Alt")
vars_variant_exonic_func$AAChange <- gsub(",$", "", vars_variant_exonic_func$AAChange, perl = TRUE)
vars_variant_exonic_func <- (vars_variant_exonic_func %>% separate_rows(AAChange, sep = ","))
vars_variant_exonic_func <- unique(vars_variant_exonic_func)

# Parse annovar output, split rows into one per transcript, and extract transcript accession IDs
vars_ann <- parse_annovar(vars_hf_hp, vars_variant_func, vars_variant_exonic_func, ensGene = ensGene, refGene = refGene, filter = c("exonic", "splicing"))


# ==================================================================================== #
# Apply somaticism filter to transcripts defined in somaticism filter transcripts file #
# ==================================================================================== #

vars_ann_som <- apply_somaticism_filter(vars_ann, somaticism_file, transcript)


# ================= #
# Apply CHIP filter #
# ================= #

# Filter for CHIP transcripts
vars_chip <- vars_ann_som[vars_ann_som[[transcript]] %in% chip_vars[[refseq_ensembl_chip_accession_column]],]

# Merge filtered variants and CHIP gene annotation data frames
vars_chip <- merge(vars_chip, chip_vars, by.x = transcript, by.y = refseq_ensembl_chip_accession_column)


# ============================================================================== #
# Filter rows if chip mutation definitions match AAChange/GeneDetail column info #
# ============================================================================== #

vars_chip_filtered <- match_mut_def(vars_chip, refseq_ensembl_suffix = refseq_ensembl_suffix)


# ==================================================== #
# Apply stricter hard filter to putative CHIP variants #
# ==================================================== #

vars_chip_filtered_put <- apply_putative_filter(vars_chip_filtered)


# ============= #
# Write to file #
# ============= #

write.csv(vars, paste(sample_id, ".all_variants.csv", sep = ""), row.names = FALSE)
write.csv(vars_hf_hp, paste(sample_id, ".all_variants.pre_filtered.csv", sep = ""), row.names = FALSE)
write.csv(vars_ann_som, paste(sample_id, ".exonic_splicing_variants.csv", sep = ""), row.names = FALSE)
write.csv(vars_chip, paste(sample_id, ".chip_transcript_variants.csv", sep = ""), row.names = FALSE)
write.csv(vars_chip_filtered, paste(sample_id, ".chip_transcript_variants.filtered.csv", sep = ""), row.names = FALSE)
write.csv(vars_chip_filtered_put, paste(sample_id, ".chip_transcript_variants.filtered.putative_filter.csv", sep = ""), row.names = FALSE)
