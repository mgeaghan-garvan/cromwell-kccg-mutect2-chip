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
    help = "Path to input VCF"
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
    c("-s", "--sample"),
    type = "character",
    default = NULL,
    help = "Sample ID"
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
    c("-c", "--chip_definitions"),
    type = "character",
    default = NULL,
    help = "Path to file containing chip definitions"
  ),
  make_option(
    c("-r", "--fasta"),
    type = "character",
    default = NULL,
    help = "Path to FASTA file"
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
  if (is.null(args$annovar) || !file.exists(args$annovar)) {
    stop("No ANNOVAR file specified")
  }
  if (is.null(args$annovar_function) || !file.exists(args$annovar_function)) {
    stop("No ANNOVAR function file specified")
  }
  if (is.null(args$annovar_exonic_function) || !file.exists(args$annovar_exonic_function)) {
    stop("No ANNOVAR exonic function file specified")
  }
  if (is.null(args$sample)) {
    stop("No sample ID specified")
  }
  if (args$gnomad_genome && args$gnomad_exome) {
    stop("Cannot specify both --gnomad_genome and --gnomad_exome")
  }
  if (!(args$gnomad_population %in% c("AF", "AF_afr", "AF_sas", "AF_amr", "AF_eas", "AF_nfe", "AF_fin", "AF_asj"))) {
    stop("Invalid gnomAD population specified")
  }
  if (is.null(args$chip_definitions) || !file.exists(args$chip_definitions)) {
    stop("No chip definitions file specified")
  }
  if (is.null(args$fasta) || !file.exists(args$fasta)) {
    stop("No FASTA file specified")
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
      "--input", "test_data/test_sample.vcf.gz",
      "--annovar", "test_data/annovar/test_sample.hg38_multianno.txt",
      "--annovar_function", "test_data/annovar/test_sample.refGene.variant_function",
      "--annovar_exonic_function", "test_data/annovar/test_sample.refGene.exonic_variant_function",
      "--sample", "test_sample",
      "--chip_definitions", "chip_mutations/chip_mutations.csv",
      "--fasta", "test_data/hg38.fa",
      "--somaticism_transcripts", "chip_mutations/somaticism_filter_transcripts.txt",
      "--output", "test_data/test_sample.chip.vcf"
    )
  )
} else {
  args <- parse_args(OptionParser(option_list = option_list))
}
check_args(args)

# ===== Load data =====

# --- Load FASTA ---
fa <- ape::read.dna(args$fasta, format = "fasta")

# --- Load VCF ---
vcf <- read.vcfR(args$input)

# --- Load CHIP definitions ---
chip_definitions <- read_csv(args$chip_definitions) %>%
  # Remove version number from transcript IDs
  mutate(
    refseq_accession = str_replace(refseq_accession, "\\.\\d+$", ""),
    ensembl_accession = str_replace(ensembl_accession, "\\.\\d+$", "")
  )
