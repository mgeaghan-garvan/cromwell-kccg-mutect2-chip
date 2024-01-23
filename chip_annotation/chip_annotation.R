#!/usr/bin/env Rscript

options(tidyverse.quiet = TRUE)

library(tidyverse, quietly = TRUE, warn.conflicts = FALSE)
library(optparse, quietly = TRUE, warn.conflicts = FALSE)

source("match_mutation.R")
# library(vcfR, quietly = TRUE, warn.conflicts = FALSE)

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
    help = "Specify that only gnomAD genome allele frequencies are present in the ANNOVAR output file; default behavior is to expect and use both gnomAD genome and exome allele frequencies, preferring exome allele frequencies when both are available; mutually exclusive with --gnomad_exome"
  ),
  make_option(
    c("-E", "--gnomad_exome"),
    action = "store_true",
    default = FALSE,
    help = "Specify that only gnomAD exome allele frequencies are present in the ANNOVAR output file; default behavior is to expect and use both gnomAD genome and exome allele frequencies, preferring exome allele frequencies when both are available; mutually exclusive with --gnomad_genome"
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
setwd("chip_annotation")
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

# ===== Function definitions =====

get_homopolymer_regions_in_sequence <- function(seq, min_size = 5) {
  # Sliding window algorithm to pick out all homopolymer regions within a sequence
  # Returns a data frame with columns "start", "end", "length", and "base"
  # Returns NULL if no homopolymer regions are found
  hpr_list <- list()
  hidx <- 1
  l <- nchar(seq)
  i <- 1
  while (i <= (l - min_size + 1)) {
    b1 <- substr(seq, i, i)
    j <- i + 1
    while (j <= l) {
      b2 <- substr(seq, j, j)
      if (b2 == b1) {
        j <- j + 1
      } else {
        break
      }
    }
    subseq_len <- j - i
    if (subseq_len >= min_size) {
      hpr_list[[hidx]] <- data.frame(start = i, end = (j - 1), length = (j - i), base = b1)
      hidx <- hidx + 1
    }
    i <- j
  }
  return(do.call(rbind, hpr_list))
}

is_altered_homopolymer_region <- function(ref_seq, alt_seq, var_pos, ref_len, alt_len) {
  # Given a reference sequence and an alternate sequence,
  # determine if the variant is within a homopolymer region
  # or otherwise alters a homopolymer region
  hp_ref <- get_homopolymer_regions_in_sequence(ref_seq)
  hp_alt <- get_homopolymer_regions_in_sequence(alt_seq)
  if (is.null(hp_ref) && is.null(hp_alt)) {
    # No homopolymer regions
    return(FALSE)
  } else if (is.null(hp_ref) && !is.null(hp_alt)) {
    # Alternate allele creates homopolymer
    return(TRUE)
  } else if (!is.null(hp_ref) && is.null(hp_alt)) {
    # Alternate allele destroys homopolymer
    return(TRUE)
  } else if (!identical(hp_ref[c("length", "base")], hp_alt[c("length", "base")])) {
    # Reference and alternate sequences have differing homopolymers (i.e. variant creates/destroys/alters a homopolymer sequence)
    return(TRUE)
  } else if (all(hp_ref$end < var_pos) && all(hp_alt$end < var_pos)) {
    # All homopolymers reside within the upstream flanking sequence (outside of variant range which starts at var_pos)
    return(FALSE)
  } else if(all(hp_ref$start > (var_pos + ref_len - 1)) && all(hp_alt$start > (var_pos + alt_len - 1))) {
    # All homopolymers reside within the downstream flanking sequence (outside of variant range which ends at (var_pos + [ref|alt]_len - 1))
    return(FALSE)
  } else {
    # All remaining variants should be within a homopolymer region
    return(TRUE)
  }
}

# ===== Load data =====

# --- Load CHIP definitions ---
chip_definitions <- read_csv(args$chip_definitions) %>%
  # Remove version number from transcript IDs
  mutate(
    refseq_accession = str_replace(refseq_accession, "\\.\\d+$", ""),
    ensembl_accession = str_replace(ensembl_accession, "\\.\\d+$", "")
  ) %>%
  # Ensure that all position columns are integers
  mutate(
    gene_genomic_start = as.integer(gene_genomic_start),
    gene_genomic_end = as.integer(gene_genomic_end),
    c_term_genomic_start = as.integer(c_term_genomic_start),
    c_term_genomic_end = as.integer(c_term_genomic_end),
  )

# --- Load sequence contexts ---
seq_contexts <- read_tsv(args$seq, col_names = FALSE)
colnames(seq_contexts) <- c("VARIANT", "SEQ")
# Values in VARIANT column have the format "chr:pos:ref:alt"
# We want to split this into separate columns
seq_contexts <- seq_contexts %>%
  separate(VARIANT, c("CHROM", "POS", "REF", "ALT"), sep = ":", remove = FALSE)
# Convert POS column to integer
seq_contexts <- seq_contexts %>%
  mutate(
    POS = as.integer(POS)
  )

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
gnomad_col_basenames <- c(
  "AF",
  "AF_popmax",
  "AF_male",
  "AF_female",
  "AF_raw",
  "AF_afr",
  "AF_sas",
  "AF_amr",
  "AF_eas",
  "AF_nfe",
  "AF_fin",
  "AF_asj",
  "AF_oth",
  "non_topmed_AF_popmax",
  "non_neuro_AF_popmax",
  "non_cancer_AF_popmax",
  "controls_AF_popmax"
)
if (args$gnomad_genome || args$gnomad_exome) {
  if (args$gnomad_genome) {
    gnomad_col_prefix <- "gnomAD_genome_"
  } else {
    gnomad_col_prefix <- "gnomAD_exome_"
  }
  gnomad_cols <- setNames(gnomad_col_basenames, paste0(gnomad_col_prefix, gnomad_cols))
} else {
  # Both genome and exome frequencies are present
  # NOTE: Should we just make this mandatory?
  if (args$gnomad_exome_first) {
    gnomad_col_prefix_1 <- "gnomAD_exome_"
    gnomad_col_prefix_2 <- "gnomAD_genome_"
  } else {
    gnomad_col_prefix_1 <- "gnomAD_genome_"
    gnomad_col_prefix_2 <- "gnomAD_exome_"
  }
  gnomad_col_regex <- paste0("^(", paste0(gnomad_col_basenames, collapse = "|"), ")\\.\\.\\.\\d+$")
  gnomad_cols <- setNames(
    grep(gnomad_col_regex, colnames(annovar), value = TRUE, perl = TRUE),
    c(
      paste0(gnomad_col_prefix_1, gnomad_col_basenames),
      paste0(gnomad_col_prefix_2, gnomad_col_basenames)
    )
  )
}

annovar <- annovar %>%
  rename(all_of(otherinfo_cols)) %>%
  rename(all_of(gnomad_cols)) %>%
  # Ensure that POS, Start, and End columns are integers
  mutate(
    POS = as.integer(POS),
    Start = as.integer(Start),
    End = as.integer(End)
  )

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
  distinct() %>%
  # Ensure that Start and End columns are integers
  mutate(
    Start = as.integer(Start),
    End = as.integer(End)
  )

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
  distinct() %>%
  # Ensure that Start and End columns are integers
  mutate(
    Start = as.integer(Start),
    End = as.integer(End)
  )

# --- Load somaticism transcripts ---
somaticism_transcripts <- read_lines(args$somaticism_transcripts)

# --- Calculate gnomAD allele frequencies ---
gnomad_genome_af_col <- paste0("gnomAD_genome_", args$gnomad_population)
gnomad_exome_af_col <- paste0("gnomAD_exome_", args$gnomad_population)
missing_vals <- c(".", "", NA)

# Create temporary gnomAD_AF_g and gnomAD_AF_e columns
annovar <- annovar %>%
  mutate(
    gnomAD_AF_g = case_when(
      !args$gnomad_exome ~ annovar[[gnomad_genome_af_col]],
      args$gnomad_exome ~ NA
    ),
    gnomAD_AF_e = case_when(
      !args$gnomad_genome ~ annovar[[gnomad_exome_af_col]],
      args$gnomad_genome ~ NA
    )
  ) %>%
  # Convert missing values to NA
  mutate(
    gnomAD_AF_g = if_else(
      gnomAD_AF_g %in% missing_vals,
      NA,
      gnomAD_AF_g
    ),
    gnomAD_AF_e = if_else(
      gnomAD_AF_e %in% missing_vals,
      NA,
      gnomAD_AF_e
    )
  ) %>%
  # Make sure gnomAD_AF_* columns are numeric
  mutate(
    gnomAD_AF_g = as.numeric(gnomAD_AF_g),
    gnomAD_AF_e = as.numeric(gnomAD_AF_e)
  ) %>%
  # Calculate gnomAD_AF column - use exome column if not NA, otherwise use genome column
  mutate(
    gnomAD_AF = if_else(
      !is.na(gnomAD_AF_e),
      gnomAD_AF_e,
      gnomAD_AF_g
    )
  ) %>%
  # Convert NAs to 0 if requested
  mutate(
    gnomAD_AF = if_else(
      is.na(gnomAD_AF),
      if_else(
        args$discard_missing_gnomad,
        NA,
        0
      ),
      gnomAD_AF
    )
  ) %>%
  # Remove temporary gnomAD_AF_g and gnomAD_AF_e columns
  select(-gnomAD_AF_g, -gnomAD_AF_e)

# --- Apply next round of hard filters ---
# Filter on gnomAD frequency
annovar <- annovar %>%
  mutate(
    HARD_FILTER_AF = case_when(
      is.na(gnomAD_AF) ~ "gnomad_af_na",
      gnomAD_AF >= 0.01 ~ "chip_gnomad_af_filter_fail",
      .default = ""
    )
  )

# Extract AD, DP, and AF (Mutect2 VAF) information from FORMAT + SAMPLE columns
# We need to use the FORMAT column
# to determine the positions of the AD, DP, and AF fields
# within the SAMPLE column, which contains semi-colon separated fields
# The AD field contains comma-separated values, one for each allele
# Similarly, the AF field contains comma-separated values, one for each alternate allele
# The DP field is a single value
annovar <- annovar %>%
  mutate(
    FORMAT_split = str_split(FORMAT, ":"),
    SAMPLE_split = str_split(SAMPLE, ":")
  ) %>%
  mutate(
    AD_i = map(FORMAT_split, ~ grep("^AD$", .x, perl = TRUE)),
    DP_i = map(FORMAT_split, ~ grep("^DP$", .x, perl = TRUE)),
    AF_i = map(FORMAT_split, ~ grep("^AF$", .x, perl = TRUE))
  ) %>%
  mutate(
    AD_REF_ALT = map2_chr(SAMPLE_split, AD_i, ~ .x[.y]),
    DP = map2_chr(SAMPLE_split, DP_i, ~ .x[.y]),
    AF = map2_chr(SAMPLE_split, AF_i, ~ .x[.y])
  ) %>%
  mutate(
    AD_REF_ALT = str_split(AD_REF_ALT, ",")
  ) %>%
  mutate(
    AD_REF = map_chr(AD_REF_ALT, ~ .x[[1]]),
    AD_ALT = map_chr(AD_REF_ALT, ~ .x[[2]])
  ) %>%
  mutate(
    AD_REF = as.integer(AD_REF),
    AD_ALT = as.integer(AD_ALT),
    DP = as.integer(DP),
    AF = as.numeric(AF)
  ) %>%
  select(-FORMAT_split, -SAMPLE_split, -AD_i, -DP_i, -AF_i, -AD_REF_ALT) %>%
  mutate(
    VAF = AD_ALT / DP
  ) %>%
  # Filter on VAF (AD/DP)
  mutate(
    HARD_FILTER_VAF = case_when(
      VAF < 0.02 ~ "chip_vaf_filter_fail",
      .default = ""
    )
  )

# --- Apply homopolymer INDEL filter ---

# Add sequence context columns to ANNOVAR table
annovar <- annovar %>%
  left_join(seq_contexts, by = c("CHROM", "POS", "REF", "ALT"))

# Ensure that the sequence context is centred on the variant and includes at least 10 bases on either side
n_additional_seq_bases <- nchar(annovar$SEQ) - nchar(annovar$REF)
uniq_n_additional_seq_bases <- unique(n_additional_seq_bases)
if (length(uniq_n_additional_seq_bases) > 1) {
  stop("Sequence context length is not consistent across variants")
}
if ((uniq_n_additional_seq_bases %% 2) != 0) {
  stop("Sequence context length (upstream + downstream) is not even")
}
# Ensure that the REF allele sequence is centred in the sequence context
seq_context_ref_start <- (uniq_n_additional_seq_bases / 2) + 1
seq_context_ref_end <- seq_context_ref_start + nchar(annovar$REF) - 1
seq_context_ref_seq <- substr(annovar$SEQ, seq_context_ref_start, seq_context_ref_end)
if (any(seq_context_ref_seq != annovar$REF)) {
  stop("REF allele sequence is not centred in the sequence context")
}
# Rename SEQ column to SEQ_REF
annovar <- annovar %>%
  rename(SEQ_REF = SEQ)
# Add SEQ_ALT column
annovar <- annovar %>%
  mutate(
    SEQ_ALT = paste(
      substr(SEQ_REF, 1, seq_context_ref_start - 1),
      ALT,
      substr(SEQ_REF, seq_context_ref_end + 1, nchar(SEQ_REF)),
      sep = ""
    )
  )
# Add HOMOPOLYMER_FILTER column
annovar <- annovar %>%
  rowwise() %>%
  mutate(
    HOMOPOLYMER_FILTER = if_else(
      is_altered_homopolymer_region(SEQ_REF, SEQ_ALT, seq_context_ref_start, nchar(REF), nchar(ALT)),
      "homopolymer_variant",
      ""
    )
  )

# --- Filter annovar annotations ---
if (args$ensembl) {
  annovar <- annovar %>% mutate(
    GeneDetail = GeneDetail.ensGene,
    AAChange = AAChange.ensGene,
    Func = Func.ensGene,
    ExonicFunc = ExonicFunc.ensGene
  )
  tx_rgx <- "(ENST\\d+(\\.\\d+)?)(:.*)?$"
} else {
  annovar <- annovar %>% mutate(
    GeneDetail = GeneDetail.refGene,
    AAChange = AAChange.refGene,
    Func = Func.refGene,
    ExonicFunc = ExonicFunc.refGene
  )
  tx_rgx <- "(NM_\\d+(\\.\\d+)?)(:.*)?$"
}

# Add a variant ID column
annovar <- annovar %>%
  mutate(VAR_ID = paste(Chr, Start, End, Ref, Alt, sep = ":"))

annovar_unfiltered <- annovar

annovar <- annovar %>%
  # Split Func column into multiple rows on semicolons
  separate_longer_delim(Func, ";") %>%
  # Filter Func for exonic and splicing variants only
  filter(
    Func %in% c("exonic", "splicing")
  ) %>%
  # Split GeneDetail and AAChange into multiple rows, one per transcript
  separate_longer_delim(GeneDetail, ";") %>%
  separate_longer_delim(AAChange, ",") %>%
  separate_longer_delim(AAChange, ";") %>%
  # Drop potential duplicates
  distinct() %>%
  # Extract the transcript ID from the GeneDetail and AAChange columns
  mutate(
    tx_gd = str_extract(GeneDetail, tx_rgx, group = 1) %>% if_else(is.na(.), "", .),
    tx_aac = str_extract(AAChange, tx_rgx, group = 1) %>% if_else(is.na(.), "", .)
  ) %>%
  mutate(
    Transcript = paste(tx_gd, tx_aac, sep = ";")
  ) %>%
  separate_longer_delim(Transcript, ";") %>%
  select(-tx_gd, -tx_aac) %>%
  filter(
    Transcript != ""
  ) %>%
  # Drop potential duplicates
  distinct() %>%
  # Add Transcript_no_version column
  mutate(
    Transcript_no_version = gsub("\\.\\d+$", "", Transcript, perl = TRUE)
  ) %>%
  # Filter for rows where AAChange and GeneDetail contains the matching transcript ID from that row
  mutate(
    Transcript_regex = paste(Transcript_no_version, "(\\.\\d+)?(:.*)?$", sep = "")
  ) %>%
  filter(
    str_detect(AAChange, Transcript_regex) | str_detect(GeneDetail, Transcript_regex)
  ) %>%
  select(-Transcript_regex) %>%
  # Drop potential duplicates
  distinct() %>%
  inner_join(annovar_function) %>%
  distinct() %>%
  # Split ExonicFunc column into multiple rows on semicolons
  separate_longer_delim(ExonicFunc, ";")

annovar <- annovar %>%
  # Join exonic variants with annovar_exonic_function table, dropping mismatches and duplicates
  filter(ExonicFunc %in% annovar_exonic_function$ExonicFunc) %>%
  left_join(annovar_exonic_function) %>%
  # Add non-exonic variants back in
  bind_rows(annovar %>% filter(!(ExonicFunc %in% annovar_exonic_function$ExonicFunc))) %>%
  distinct()

# Annotate variants that were filtered out
annovar_filtered <- annovar_unfiltered %>%
  anti_join(annovar, by = "VAR_ID") %>%
  mutate(
    EXONIC_SPLICING_FILTER = "not_exonic_or_splicing_variant"
  )

# --- Apply somaticism filter ---
annovar <- annovar %>%
  mutate(ad_dp_binom_test = map2_dbl(AD_ALT, DP, ~ binom.test(.x, .y, 0.5)$p.value)) %>%
  mutate(
    SOMATICISM_FILTER = case_when(
      Transcript_no_version %in% somaticism_transcripts & ad_dp_binom_test >= 0.001 ~ "chip_somaticism_filter_fail",
      .default = ""
    )
  )

# --- Filter for CHIP transcripts ---
annovar_chip_tx_unfiltered <- annovar

if (args$ensembl) {
  annovar <- annovar %>%
    filter(
      Transcript_no_version %in% chip_definitions$ensembl_accession
    )
} else {
  annovar <- annovar %>%
    filter(
      Transcript_no_version %in% chip_definitions$refseq_accession
    )
}

# Annotate variants that were filtered out
annovar_chip_tx_filtered <- annovar_chip_tx_unfiltered %>%
  anti_join(annovar, by = "VAR_ID") %>%
  mutate(
    CHIP_TRANSCRIPT_FILTER = "exonic_or_splicing_variant_not_in_chip_transcript"
  )

# Join ANNOVAR table with CHIP definitions
if (args$ensembl) {
  annovar <- annovar %>%
    left_join(chip_definitions, by = c("Transcript_no_version" = "ensembl_accession"), relationship = "many-to-many")
} else {
  annovar <- annovar %>%
    left_join(chip_definitions, by = c("Transcript_no_version" = "refseq_accession"), relationship = "many-to-many")
}

# Drop potential duplicates
annovar <- annovar %>%
  distinct()

# --- Match mutations to CHIP definitions ---
annovar_tmp <- annovar %>%
  mutate(
    mut_in_c_term = case_when(
      (
        Chr == chr & (
          (c_term_genomic_start <= Start & Start <= c_term_genomic_end) |
          (c_term_genomic_start <= End & End <= c_term_genomic_end) |
          (Start <= c_term_genomic_start & c_term_genomic_end <= End)
        )
      ) ~ TRUE,
      .default = FALSE
    )
  ) %>%
  rowwise() %>%
  mutate(
    gene_detail_info = list(parse_gene_detail(GeneDetail)),
    aa_change_info = list(parse_aa_change(AAChange, ExonicFunc)),
    chip_info = list(parse_chip_def(mutation_definition, mutation_class))
  ) %>%
  mutate(
    variant_class = case_when(
      Func == "splicing" ~ "splice_site",
      (
        Func != "splicing" &&
        aa_change_info$mutation_class == "nonsynonymous" &&
        mutation_class == "missense" &&
        !is.na(aa_change_info$nonsynonymous_type) &&
        aa_change_info$nonsynonymous_type == "sub" &&
        !is.na(aa_change_info$mutation_protein_info$ref_start) &&
        !is.na(aa_change_info$mutation_protein_info$alt)
      ) ~ "missense",
      .default = aa_change_info$mutation_class
    )
  ) %>%
  mutate(
    # Dynamically picks the appropriate function to use based on the mutation pattern
    chip_mutation_match_filter = chip_def_match_funcs[[chip_info$mutation_pattern]](
      chip_info = chip_info,
      aa_change_info = aa_change_info,
      gene_detail_info = gene_detail_info,
      mut_class = variant_class,
      chip_mut_c_term_start = c_term_genomic_start,
      chip_mut_c_term_end = c_term_genomic_end,
      mut_start = Start,
      mut_end = End,
      mut_in_c_term = mut_in_c_term
    )
  ) %>%
  ungroup() %>%
  mutate(
    CHIP_MUTATION_FILTER = case_when(
      is.na(variant_class) ~ "chip_mutation_match_filter_fail",
      variant_class != mutation_class ~ "chip_mutation_match_filter_fail",
      .default = chip_mutation_match_filter
    )
  ) %>%
  select(-gene_detail_info, -aa_change_info, -chip_info, -variant_class, -chip_mutation_match_filter, -mut_in_c_term)

# --- Apply a putative CHIP filter ---


# --- Finalise combined filter ---


# --- Construct output VCF for annotation ---


# --- Write output VCF ---


# --- Write output CSVs ---


# --- Load VCFs ---
# input_vcf <- read.vcfR(args$input)
# if (!is.null(args$failed_variants) && file.exists(args$failed_variants)) {
#   failed_variants_vcf <- read.vcfR(args$failed_variants)
# }
