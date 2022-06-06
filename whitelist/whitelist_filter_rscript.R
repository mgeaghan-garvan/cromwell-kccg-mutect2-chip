# library(data.table, quietly = T)
library(tidyr)
library(stringr)

# Source R scripts
source("./import/gnomad.R")
source("./import/hard_filter.R")
source("./import/homopolymers.R")
source("./import/parse_annovar.R")
source("./import/match_mutation.R")
source("./import/apply_putative_filter.R")


# ========================== #
# Get command-line arguments #
# ========================== #

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9) {
  stop(paste(
    "Incorrect number of arguments!",
    "Usage: Rscript whitelist_filter_rscript.R <ANNOVAR OUTPUT TABLE> <ANNOVAR OUTPUT VCF> <ANNOVAR VARIANT FUNCTION TABLE> <TUMOR SAMPLE NAME> <GNOMAD SOURCE> <GNOMAD SUBPOPULATION CODE> <TREAT MISSING AF AS RARE> <CHIP DEFINITION FILE> <TRANSCRIPT PROTEIN LENGTHS FILE> <FASTA REFERENCE FILE>",
    "    ANNOVAR OUTPUT TABLE/VCF:         output txt and vcf files from Annovar.",
    "    ANNOVAR VARIANT FUNCTION TABLE:   variant function output file from Annovar.",
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
annovar_variant_func_out <- args[3]
tumor_sample_name <- args[4]
gnomad_source <- args[5]
gnomad_pop <- args[6]
treat_missing_as_rare <- args[7]
chip_def_file <- args[8]
transcript_prot_file <- args[9]
fasta_file <- args[10]


# ============================ #
# Check command-line arguments #
# ============================ #

if (!file.exists(annovar_text_out)) {
  stop("Input annovar table file does not exist.")
}
if (!file.exists(annovar_vcf_out)) {
  stop("Input annovar vcf file does not exist.")
}
if (!file.exists(annovar_variant_func_out)) {
  stop("Input annovar variant function file does not exist.")
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


# =========================================== #
# Select between using RefSeq and Ensembl IDs #
# =========================================== #

ensembl_refseq <- "refseq"  # Use Ensembl IDs by default. Change this to "refseq" if you want to use RefSeq IDs
ensGene <- TRUE
refGene <- FALSE
refseq_ensembl_suffix <- "ensGene"
refseq_ensembl_chip_accession_column <- "ensembl_accession"
refseq_ensembl_variant_func_regex_sub <- "(ENST\\d+)([^\\d]|$)"
refseq_ensembl_variant_func_regex_match <- "^ENST\\d+$"
if (ensembl_refseq == "refseq") {
  ensGene <- FALSE
  refGene <- TRUE
  refseq_ensembl_suffix <- "refGene"
  refseq_ensembl_chip_accession_column <- "refseq_accession"
  refseq_ensembl_variant_func_regex_sub <- "(NM_\\d+)([^\\d]|$)"
  refseq_ensembl_variant_func_regex_match <- "^NM_\\d+$"
}

aachange <- paste("AAChange", refseq_ensembl_suffix, sep = ".")
genedetail <- paste("GeneDetail", refseq_ensembl_suffix, sep = ".")
transcript <- paste("Transcript", refseq_ensembl_suffix, sep = ".")
exonic_func <- paste("ExonicFunc", refseq_ensembl_suffix, sep = ".")


# ================================================ #
# Apply filters to variants in homopolymer regions #
# ================================================ #

# NOTE: For now, we'll just be doing INDELs (HP_SNVS = FALSE), but in the future we should make this configurable via CLI arguments
vars_hf_hp <- apply_homopolymer_filter(vars_hf, fasta_file, exonic_func_col = exonic_func, HP_SNVS = FALSE)


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

# Parse annovar output, split rows into one per transcript, and extract transcript accession IDs
vars_ann <- parse_annovar(vars_hf_hp, vars_variant_func, ensGene = ensGene, refGene = refGene, filter = c("exonic", "splicing"))

# Filter for CHIP transcripts
vars_chip <- vars_ann[vars_ann[[transcript]] %in% chip_vars[[refseq_ensembl_chip_accession_column]],]

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
write.csv(vars_ann, paste(sample_id, ".exonic_splicing_variants.csv", sep = ""), row.names = FALSE)
write.csv(vars_chip, paste(sample_id, ".chip_transcript_variants.csv", sep = ""), row.names = FALSE)
write.csv(vars_chip_filtered, paste(sample_id, ".chip_transcript_variants.filtered.csv", sep = ""), row.names = FALSE)
write.csv(vars_chip_filtered_put, paste(sample_id, ".chip_transcript_variants.filtered.putative_filter.csv", sep = ""), row.names = FALSE)
