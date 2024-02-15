#!/usr/bin/env Rscript

options(tidyverse.quiet = TRUE)

library(tidyverse, quietly = TRUE, warn.conflicts = FALSE)
library(optparse, quietly = TRUE, warn.conflicts = FALSE)

# Parse command line arguments
option_list <- list(
  make_option(
    c("-i", "--input_vcf"),
    type = "character",
    default = NULL,
    help = "Path to the input VCF file"
  ),
  make_option(
    c("-t", "--prevalence_threshold"),
    type = "numeric",
    default = NULL,
    help = "Prevalence threshold for filtering variants"
  )
)

check_args <- function(args) {
  if (is.null(args$input_vcf) || !file.exists(args$input_vcf)) {
    stop("No VCF file provided or file does not exist")
  }
  if (is.null(args$prevalence_threshold)) {
    stop("Please provide the prevalence threshold for filtering variants")
  }
  if (!is.numeric(args$prevalence_threshold)) {
    stop("Prevalence threshold must be a number")
  }
  if (args$prevalence_threshold < 0 || args$prevalence_threshold > 1) {
    stop("Prevalence threshold must be between 0 and 1")
  }
}

args <- parse_args(OptionParser(option_list = option_list))
check_args(args)

# Read VCF file header
vcf_header <- read_lines(args$input_vcf) %>%
  keep(~ str_detect(., "^#"))
# Remove the the FORMAT and SAMPLE columns from the header
vcf_header_no_format <- vcf_header
vcf_header_no_format[length(vcf_header_no_format)] <- gsub(
  "\\tFORMAT.*$",
  "",
  vcf_header_no_format[length(vcf_header_no_format)]
)

# Get VCF header columns
vcf_header_cols <- vcf_header %>%
  tail(1) %>%
  str_split("\t") %>%
  unlist() %>%
  gsub("^#", "", .)

# Read VCF file body
vcf <- read_tsv(
  args$input_vcf,
  comment = "#",
  col_names = vcf_header_cols
)

# Filter variants based on prevalence
vcf_filtered <- vcf %>%
  mutate(
    across(
      all_of(10:ncol(vcf)),
      ~ str_split(.x, "/|\\|")
    )
  ) %>%
  mutate(
    across(
      all_of(10:ncol(vcf)),
      ~ map_lgl(., ~ any(!(.x %in% c("0", "."))))
    )
  ) %>%
  mutate(
    prevalence = rowSums(
      across(
        all_of(10:ncol(vcf)),
        ~ .x
      )
    ) / (ncol(vcf) - 9)
  ) %>%
  mutate(
    FILTER = if_else(
      prevalence >= args$prevalence_threshold,
      "chip_cohort_prevalence_filter_fail",
      "."
    ),
    INFO = paste0("CHIP_Cohort_Prevalence=", prevalence)
  ) %>%
  # Annotate multi-allelic variants with a CHIP_Multiallelic_Cohort_Prevalence_Filters INFO field
  group_by(
    CHROM,
    POS,
    REF
  ) %>%
  mutate(
    INFO_CHIP_Multiallelic_Cohort_Prevalence_Filters = case_when(
      length(unique(FILTER)) > 1 & FILTER == "chip_cohort_prevalence_filter_fail" ~ "CHIP_Multiallelic_Cohort_Prevalence_Filters=chip_cohort_prevalence_filter_fail",
      .default = ""
    )
  ) %>%
  mutate(
    INFO = case_when(
      INFO_CHIP_Multiallelic_Cohort_Prevalence_Filters != "" ~ paste0(INFO, ";", INFO_CHIP_Multiallelic_Cohort_Prevalence_Filters),
      .default = INFO
    )
  ) %>%
  ungroup() %>%
  select(
    CHROM,
    POS,
    ID,
    REF,
    ALT,
    QUAL,
    FILTER,
    INFO
  ) %>%
  distinct()

filter_header = paste0("##FILTER=<ID=chip_cohort_prevalence_filter_fail,Description=\"Variant fails CHIP cohort prevalence filter with threshold ", prevalence, "\">")
info_header = c(
  "##INFO=<ID=CHIP_Cohort_Prevalence,Number=A,Type=Float,Description=\"Prevalence of variant in the CHIP cohort\">",
  "##INFO=<ID=CHIP_Multiallelic_Cohort_Prevalence_Filters,Number=A,Type=String,Description=\"Allele-specific CHIP cohort prevalence filter for multiallelic variants\">"
)

# Insert new FILTER and INFO fields into header
vcf_header_no_format <- vcf_header_no_format %>%
  append(
    c(filter_header, info_header),
    after = length(vcf_header_no_format) - 1
  )

# Write filtered VCF file
output_vcf <- gsub("\\.vcf(\\.gz)?$", "annotated.vcf", args$input_vcf)
vcf_header_no_format %>%
  write_lines(output_vcf)

vcf_filtered %>%
  write_tsv(output_vcf, append = TRUE, col_names = FALSE)
