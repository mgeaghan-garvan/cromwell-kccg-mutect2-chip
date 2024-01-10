#!/usr/bin/env Rscript

options(tidyverse.quiet = TRUE)

library(tidyverse, quietly = TRUE, warn.conflicts = FALSE)
library(optparse, quietly = TRUE, warn.conflicts = FALSE)
library(vcfR, quietly = TRUE, warn.conflicts = FALSE)

# Parse command line arguments
option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    default = NULL,
    help = "Path to input VCF containing variants to compare against known CHIP mutations"
  ),
  make_option(
    c("-I", "--failed_variants"),
    type = "character",
    default = NULL,
    help = "Path to an optional VCF containing variants that failed to pass filters; FILTER and INFO fields will be copied to the output VCF"
  ),
  make_option(
    c("-s", "--sample"),
    type = "character",
    default = NULL,
    help = "Sample ID"
  ),
  make_option(
    c("-c", "--chip_definitions"),
    type = "character",
    default = NULL,
    help = "Path to file containing chip definitions"
  ),
  make_option(
    c("-q", "--seq"),
    type = "character",
    default = NULL,
    help = "Path to CSV file containing sequence contexts for each variant; columns must be 'CHROM', 'POS', 'REF', 'ALT', 'SEQ'; sequence must be centered on the variant and include at least 10 bases on either side of the variant"
  ),
  make_option(
    c("-a", "--annovar"),
    type = "character",
    default = NULL,
    help = "Path to ANNOVAR table containing refGene annotations and gnomAD allele frequencies"
  ),
  make_option(
    c("-f", "--annovar_function"),
    type = "character",
    default = NULL,
    help = "Path to ANNOVAR intermediate file containing variant function annotations"
  ),
  make_option(
    c("-F", "--annovar_exonic_function"),
    type = "character",
    default = NULL,
    help = "Path to ANNOVAR intermediate file containing variant exonic function annotations"
  ),
  make_option(
    c("-e", "--ensembl"),
    action = "store_true",
    default = FALSE,
    help = "Use Ensembl transcript IDs instead of RefSeq transcript IDs"
  ),
  make_option(
    c("-G", "--gnomad_genome"),
    action = "store_true",
    default = FALSE,
    help = "Only use gnomAD genome allele frequencies; default behavior is to use both gnomAD genome and exome allele frequencies, preferring exome allele frequencies when both are available; mutually exclusive with --gnomad_exome"
  ),
  make_option(
    c("-E", "--gnomad_exome"),
    action = "store_true",
    default = FALSE,
    help = "Only use gnomAD exome allele frequencies; default behavior is to use both gnomAD genome and exome allele frequencies, preferring exome allele frequencies when both are available; mutually exclusive with --gnomad_genome"
  ),
  make_option(
    c("-x", "--gnomad_exome_first"),
    action = "store_true",
    default = FALSE,
    help = "When both gnomAD genome and exome allele frequencies are available, specify that the exome frequencies come first in the column order; by default genome frequencies are expected first; column order of gnomAD frequencies is determined by their order in the ANNOVAR command; incompatible with --gnomad_genome and --gnomad_exome"
  ),
  make_option(
    c("-M", "--discard_missing_gnomad"),
    action = "store_true",
    default = FALSE,
    help = "Discard variants with missing gnomAD allele frequencies; default behavior is to assign missing gnomAD allele frequencies a value of 0"
  ),
  make_option(
    c("-P", "--gnomad_population"),
    type = "character",
    default = "AF",
    help = "gnomAD population to use for allele frequencies; default behavior is to use all populations ('AF'); options are: 'AF', 'AF_afr', 'AF_sas', 'AF_amr', 'AF_eas', 'AF_nfe', 'AF_fin', 'AF_asj'"
  ),
  make_option(
    c("-S", "--somaticism_transcripts"),
    type = "character",
    default = NULL,
    help = "Path to file containing somaticism transcript IDs"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = NULL,
    help = "Path to output file"
  )
)

check_args <- function(args) {
  if (is.null(args$input) || !file.exists(args$input)) {
    stop("No input file specified")
  }
  if (is.null(args$sample)) {
    stop("No sample ID specified")
  }
  if (is.null(args$chip_definitions) || !file.exists(args$chip_definitions)) {
    stop("No chip definitions file specified")
  }
  if (is.null(args$seq) || !file.exists(args$seq)) {
    stop("No sequence context file specified")
  }
  if (is.null(args$annovar) || !file.exists(args$annovar)) {
    stop("No ANNOVAR file specified")
  }
  if (is.null(args$annovar_function) || !file.exists(args$annovar_function)) {
    stop("No ANNOVAR function file specified")
  }
  if (is.null(args$annovar_exonic_function) || !file.exists(args$annovar_exonic_function)) {
    stop("No ANNOVAR exonic function file specified")
  }
  if (args$gnomad_genome && args$gnomad_exome) {
    stop("Cannot specify both --gnomad_genome and --gnomad_exome")
  }
  if (args$gnomad_exome_first && (args$gnomad_genome || args$gnomad_exome)) {
    stop("Cannot specify --gnomad_exome_first with --gnomad_genome or --gnomad_exome")
  }
  if (!(args$gnomad_population %in% c("AF", "AF_afr", "AF_sas", "AF_amr", "AF_eas", "AF_nfe", "AF_fin", "AF_asj"))) {
    stop("Invalid gnomAD population specified")
  }
  if (is.null(args$somaticism_transcripts) || !file.exists(args$somaticism_transcripts)) {
    stop("No somaticism transcripts file specified")
  }
  if (is.null(args$output)) {
    stop("No output file specified")
  }
}

DEBUG <- TRUE
if (DEBUG) {
  args <- parse_args(
    OptionParser(option_list = option_list),
    c(
      "--input", "test_data/mutect2_out/test_sample.vcf.gz",
      "--failed_variants", "test_data/filter/test_sample.non_chip_genes.so.chip_filter.vcf",
      "--sample", "test_sample",
      "--chip_definitions", "test_data/chip_mutations/chip_mutations.chr.csv",
      "--seq", "test_data/filter/test_sample.chip_genes.norm.seq.tsv",
      "--annovar", "test_data/annovar/test_sample.chip_genes.norm.no_info.annot.hg38_multianno.txt",
      "--annovar_function", "test_data/annovar/test_sample.chip_genes.norm.no_info.annot.refGene.variant_function",
      "--annovar_exonic_function", "test_data/annovar/test_sample.chip_genes.norm.no_info.annot.refGene.exonic_variant_function",
      "--somaticism_transcripts", "chip_mutations/somaticism_filter_transcripts.txt",
      "--output", "test_data/test_sample.chip_annotations.vcf"
    )
  )
} else {
  args <- parse_args(OptionParser(option_list = option_list))
}
check_args(args)

# ===== Load data =====

# --- Load CHIP definitions ---
chip_definitions <- read_csv(args$chip_definitions) %>%
  # Remove version number from transcript IDs
  mutate(
    refseq_accession = str_replace(refseq_accession, "\\.\\d+$", ""),
    ensembl_accession = str_replace(ensembl_accession, "\\.\\d+$", "")
  )

# --- Load sequence contexts ---
seq_contexts <- read_tsv(args$seq, col_names = FALSE)
colnames(seq_contexts) <- c("VARIANT", "SEQ")
# Values in VARIANT column have the format "chr:pos:ref:alt"
# We want to split this into separate columns
seq_contexts <- seq_contexts %>%
  separate(VARIANT, c("CHROM", "POS", "REF", "ALT"), sep = ":", remove = FALSE)

# --- Load ANNOVAR annotations ---
annovar <- read_tsv(args$annovar)

# Get new names for "Otherinfo" columns generated by ANNOVAR
otherinfo_cols <- c(
  ANNOVAR_alt_af = "Otherinfo1",
  ANNOVAR_qual = "Otherinfo2",
  ANNOVAR_alt_ad = "Otherinfo3",
  CHROM = "Otherinfo4",
  POS = "Otherinfo5",
  ID = "Otherinfo6",
  REF = "Otherinfo7",
  ALT = "Otherinfo8",
  QUAL = "Otherinfo9",
  FILTER = "Otherinfo10",
  INFO = "Otherinfo11",
  FORMAT = "Otherinfo12",
  SAMPLE = "Otherinfo13"
)

# Determine new column names for gnomAD columns
if (args$genome) {
  gnomad_cols <- c()
} else if (args$exome) {
  gnomad_cols <- c()
} else {
  # Both genome and exome frequencies are present
  # NOTE: Should we just make this mandatory?
  gnomad_cols <- c()
}

annovar <- annovar %>%
  rename(all_of(otherinfo_cols))

# Load ANNOVAR intermediate files
# First, generate regex strings for transcript IDs
if (args$ensembl) {
  transcript_id_variant_function_regex <- "^(ENST_\\d+)([^\\d].*)?$$"
  transcript_id_regex <- "^ENST_\\d+$"
} else {
  transcript_id_variant_function_regex <- "^(NM_\\d+)([^\\d].*)?$"
  transcript_id_regex <- "^NM_\\d+$"
}

annovar_function <- read_tsv(args$annovar_function, col_names = FALSE) %>%
  # Select just the first 7 columns and remove duplicates
  select(1:7) %>%
  distinct() %>%
  # Rename columns
  rename(
    Func = X1,
    Transcript = X2,
    Chr = X3,
    Start = X4,
    End = X5,
    Ref = X6,
    Alt = X7
  ) %>%
  # Clean up Transcript column
  # Values will look like "NM_021170(dist=2318),NM_005101(dist=11082)"
  # We want to only retain the transcript IDs
  separate_longer_delim(Transcript, ",") %>%
  mutate(
    Transcript = str_replace(Transcript, transcript_id_variant_function_regex, "\\1")
  ) %>%
  filter(
    str_detect(Transcript, transcript_id_regex)
  ) %>%
  distinct()

annovar_exonic_function <- read_tsv(args$annovar_exonic_function, col_names = FALSE) %>%
  # Select columns 2-8 and remove duplicates
  select(2:8) %>%
  distinct() %>%
  # Rename columns
  rename(
    ExonicFunc = X2,
    AAChange = X3,
    Chr = X4,
    Start = X5,
    End = X6,
    Ref = X7,
    Alt = X8
  ) %>%
  # Clean up AAChange column by removing trailing commas
  mutate(
    AAChange = str_replace(AAChange, ",$", "")
  ) %>%
  # Split AAChange column into multiple rows, split on commas
  separate_longer_delim(AAChange, ",") %>%
  distinct()

# --- Load somaticism transcripts ---
somaticism_transcripts <- read_lines(args$somaticism_transcripts)

# --- Load VCF ---
# vcf <- read.vcfR(args$input)
