#!/usr/bin/env Rscript

options(tidyverse.quiet = TRUE)

library(tidyverse, quietly = TRUE, warn.conflicts = FALSE)
library(optparse, quietly = TRUE, warn.conflicts = FALSE)

source("/R/scripts/match_mutation.R")

# Parse command line arguments
option_list <- list(
  make_option(
    c("-v", "--vcf_header"),
    type = "character",
    default = NULL,
    help = "Path to the header of the VCF being annotated"
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
    c("-o", "--output_prefix"),
    type = "character",
    default = NULL,
    help = "Prefix for output files"
  )
)

check_args <- function(args) {
  if (is.null(args$vcf_header) || !file.exists(args$vcf_header)) {
    stop("No VCF header file specified")
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
  if (is.null(args$output_prefix)) {
    stop("No output prefix specified")
  }
}

args <- parse_args(OptionParser(option_list = option_list))
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
  } else if (
    all(hp_ref$end < var_pos | hp_ref$start > (var_pos + ref_len - 1)) &&
      all(hp_alt$end < var_pos | hp_alt$start > (var_pos + alt_len - 1))
  ) {
    # All homopolymers reside within the upstream or downstream flanking sequences
    # (outside of variant range which starts at var_pos and ends at (var_pos + [ref|alt]_len - 1)
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

# --- Create main data frame ---
df <- annovar

# --- Calculate gnomAD allele frequencies ---
gnomad_genome_af_col <- paste0("gnomAD_genome_", args$gnomad_population)
gnomad_exome_af_col <- paste0("gnomAD_exome_", args$gnomad_population)
missing_vals <- c(".", "", NA)

# Create temporary gnomAD_AF_g and gnomAD_AF_e columns
df <- df %>%
  mutate(
    gnomAD_AF_g = case_when(
      !args$gnomad_exome ~ df[[gnomad_genome_af_col]],
      args$gnomad_exome ~ NA
    ),
    gnomAD_AF_e = case_when(
      !args$gnomad_genome ~ df[[gnomad_exome_af_col]],
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
df <- df %>%
  mutate(
    HARD_FILTER_AF = case_when(
      is.na(gnomAD_AF) ~ "gnomad_af_na",
      gnomAD_AF >= 0.001 ~ "chip_gnomad_af_filter_fail",
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
df <- df %>%
  mutate(
    FORMAT_split = str_split(FORMAT, ":"),
    SAMPLE_split = str_split(SAMPLE, ":")
  ) %>%
  mutate(
    AD_i = map(FORMAT_split, ~ grep("^AD$", .x, perl = TRUE)),
    DP_i = map(FORMAT_split, ~ grep("^DP$", .x, perl = TRUE)),
    AF_i = map(FORMAT_split, ~ grep("^AF$", .x, perl = TRUE)),
    F1R2_i = map(FORMAT_split, ~ grep("^F1R2$", .x, perl = TRUE)),
    F2R1_i = map(FORMAT_split, ~ grep("^F2R1$", .x, perl = TRUE))
  ) %>%
  mutate(
    AD_REF_ALT = map2_chr(SAMPLE_split, AD_i, ~ .x[.y]),
    DP = map2_chr(SAMPLE_split, DP_i, ~ .x[.y]),
    AF = map2_chr(SAMPLE_split, AF_i, ~ .x[.y]),
    F1R2_REF_ALT = map2_chr(SAMPLE_split, F1R2_i, ~ .x[.y]),
    F2R1_REF_ALT = map2_chr(SAMPLE_split, F2R1_i, ~ .x[.y])
  ) %>%
  mutate(
    AD_REF_ALT = str_split(AD_REF_ALT, ","),
    F1R2_REF_ALT = str_split(F1R2_REF_ALT, ","),
    F2R1_REF_ALT = str_split(F2R1_REF_ALT, ",")
  ) %>%
  mutate(
    AD_REF = map_chr(AD_REF_ALT, ~ .x[[1]]),
    AD_ALT = map_chr(AD_REF_ALT, ~ .x[[2]]),
    F1R2_REF = map_chr(F1R2_REF_ALT, ~ .x[[1]]),
    F1R2_ALT = map_chr(F1R2_REF_ALT, ~ .x[[2]]),
    F2R1_REF = map_chr(F2R1_REF_ALT, ~ .x[[1]]),
    F2R1_ALT = map_chr(F2R1_REF_ALT, ~ .x[[2]])
  ) %>%
  mutate(
    AD_REF = as.integer(AD_REF),
    AD_ALT = as.integer(AD_ALT),
    DP = as.integer(DP),
    AF = as.numeric(AF),
    F1R2_REF = as.integer(F1R2_REF),
    F1R2_ALT = as.integer(F1R2_ALT),
    F2R1_REF = as.integer(F2R1_REF),
    F2R1_ALT = as.integer(F2R1_ALT)
  ) %>%
  select(
    -FORMAT_split,
    -SAMPLE_split,
    -AD_i,
    -DP_i,
    -AF_i,
    -AD_REF_ALT,
    -F1R2_i,
    -F2R1_i,
    -F1R2_REF_ALT,
    -F2R1_REF_ALT
  ) %>%
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
df <- df %>%
  left_join(seq_contexts, by = c("CHROM", "POS", "REF", "ALT"))

# Ensure that the sequence context is centred on the variant and includes at least 10 bases on either side
n_additional_seq_bases <- nchar(df$SEQ) - nchar(df$REF)
uniq_n_additional_seq_bases <- unique(n_additional_seq_bases)
if (length(uniq_n_additional_seq_bases) > 1) {
  stop("Sequence context length is not consistent across variants")
}
if ((uniq_n_additional_seq_bases %% 2) != 0) {
  stop("Sequence context length (upstream + downstream) is not even")
}
# Ensure that the REF allele sequence is centred in the sequence context
seq_context_ref_start <- (uniq_n_additional_seq_bases / 2) + 1
seq_context_ref_end <- seq_context_ref_start + nchar(df$REF) - 1
seq_context_ref_seq <- substr(df$SEQ, seq_context_ref_start, seq_context_ref_end)
if (any(seq_context_ref_seq != df$REF)) {
  stop("REF allele sequence is not centred in the sequence context")
}
# Rename SEQ column to SEQ_REF
df <- df %>%
  rename(SEQ_REF = SEQ)
# Add SEQ_ALT column
df <- df %>%
  mutate(
    SEQ_ALT = paste(
      substr(SEQ_REF, 1, seq_context_ref_start - 1),
      ALT,
      substr(SEQ_REF, seq_context_ref_end + 1, nchar(SEQ_REF)),
      sep = ""
    )
  )
# Add HOMOPOLYMER_FILTER column
df <- df %>%
  rowwise() %>%
  mutate(
    HOMOPOLYMER_FILTER = if_else(
      is_altered_homopolymer_region(SEQ_REF, SEQ_ALT, seq_context_ref_start, nchar(REF), nchar(ALT)),
      "homopolymer_variant",
      ""
    )
  ) %>%
  # Add HOMOPOLYMER_INDEL_FILTER column
  mutate(
    HOMOPOLYMER_INDEL_FILTER = case_when(
      (
        HOMOPOLYMER_FILTER == "homopolymer_variant" &
          (nchar(REF) > 1 | nchar(ALT) > 1) &
          (AD_ALT < 10 | VAF < 0.1)
      ) ~ "chip_homopolymer_indel_filter_fail",
      .default = ""
    )
  )

# --- Filter annovar annotations ---
if (args$ensembl) {
  df <- df %>% mutate(
    GeneDetail = GeneDetail.ensGene,
    AAChange = AAChange.ensGene,
    Func = Func.ensGene,
    ExonicFunc = ExonicFunc.ensGene
  )
  tx_rgx <- "(ENST\\d+(\\.\\d+)?)(:.*)?$"
} else {
  df <- df %>% mutate(
    GeneDetail = GeneDetail.refGene,
    AAChange = AAChange.refGene,
    Func = Func.refGene,
    ExonicFunc = ExonicFunc.refGene
  )
  tx_rgx <- "(NM_\\d+(\\.\\d+)?)(:.*)?$"
}

# Add a variant ID column
df <- df %>%
  mutate(VAR_ID = paste(Chr, Start, End, Ref, Alt, sep = ":"))

df_unfiltered <- df

df <- df %>%
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

df <- df %>%
  # Join exonic variants with annovar_exonic_function table, dropping mismatches and duplicates
  filter(ExonicFunc %in% annovar_exonic_function$ExonicFunc) %>%
  left_join(annovar_exonic_function) %>%
  # Add non-exonic variants back in
  bind_rows(df %>% filter(!(ExonicFunc %in% annovar_exonic_function$ExonicFunc))) %>%
  distinct()

# Annotate variants that were filtered out
df_filtered <- df_unfiltered %>%
  anti_join(df, by = "VAR_ID") %>%
  mutate(
    EXONIC_SPLICING_FILTER = "not_exonic_or_splicing_variant"
  )

# --- Apply somaticism filter ---
df <- df %>%
  mutate(
    ad_dp_binom_test = map2_dbl(AD_ALT, DP, ~ if (.y > 0) { binom.test(.x, .y, 0.5)$p.value } else { 1 })
  ) %>%
  mutate(
    SOMATICISM_FILTER = case_when(
      Transcript_no_version %in% somaticism_transcripts & ad_dp_binom_test >= 0.001 ~ "chip_somaticism_filter_fail",
      .default = ""
    )
  )

# --- Filter for CHIP transcripts ---
df_chip_tx_unfiltered <- df

if (args$ensembl) {
  df <- df %>%
    filter(
      Transcript_no_version %in% chip_definitions$ensembl_accession
    )
} else {
  df <- df %>%
    filter(
      Transcript_no_version %in% chip_definitions$refseq_accession
    )
}

# Annotate variants that were filtered out
df_chip_tx_filtered <- df_chip_tx_unfiltered %>%
  anti_join(df, by = "VAR_ID") %>%
  mutate(
    CHIP_TRANSCRIPT_FILTER = "exonic_or_splicing_variant_not_in_chip_transcript"
  )

# Join ANNOVAR table with CHIP definitions
if (args$ensembl) {
  df <- df %>%
    left_join(chip_definitions, by = c("Transcript_no_version" = "ensembl_accession"), relationship = "many-to-many")
} else {
  df <- df %>%
    left_join(chip_definitions, by = c("Transcript_no_version" = "refseq_accession"), relationship = "many-to-many")
}

# Drop potential duplicates
df <- df %>%
  distinct()

# --- Match mutations to CHIP definitions ---
df <- df %>%
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
    variant_class = case_when(
      is.na(variant_class) ~ "",
      .default = variant_class
    )
  ) %>%
  mutate(
    # Dynamically picks the appropriate function to use based on the mutation pattern
    CHIP_MUTATION_FILTER = case_when(
      is.na(variant_class) ~ "chip_mutation_match_filter_fail",
      variant_class == "" ~ "chip_mutation_match_filter_fail",
      variant_class != mutation_class ~ "chip_mutation_match_filter_fail",
      .default = chip_def_match_funcs[[chip_info$mutation_pattern]](
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
    )
  ) %>%
  ungroup() %>%
  select(-gene_detail_info, -aa_change_info, -chip_info, -variant_class, -mut_in_c_term)

# --- Apply a putative CHIP filter ---
df <- df %>%
  mutate(
    PUTATIVE_CHIP_FILTER = case_when(
      !is.na(putative) &
        putative == TRUE &
        (
          AD_ALT < 5 |
          VAF > 0.2 |
          F1R2_REF < 2 |
          F1R2_ALT < 2 |
          F2R1_REF < 2 |
          F2R1_ALT < 2
        ) ~ "chip_putative_filter_fail",
      .default = ""
    )
  )

# --- Finalise combined filter and CHIP info ---
df <- df %>%
  # Combine variant-level filters into a single semicolon-separated string
  group_by(
    CHROM,
    POS,
    REF,
    ALT
  ) %>%
  mutate(
    VARIANT_LEVEL_FILTER_STR_COMB = paste(
      FILTER,
      HARD_FILTER_AF,
      HARD_FILTER_VAF,
      HOMOPOLYMER_INDEL_FILTER,
      SOMATICISM_FILTER,
      sep = ";"
    )
  ) %>%
  mutate(
    VARIANT_LEVEL_FILTER_STR_COMB = paste(VARIANT_LEVEL_FILTER_STR_COMB, collapse = ";")
  ) %>%
  ungroup() %>%
  mutate(
    VARIANT_LEVEL_FILTER_LIST = str_split(VARIANT_LEVEL_FILTER_STR_COMB, ";")
  ) %>%
  mutate(
    VARIANT_LEVEL_FILTER_LIST = map(VARIANT_LEVEL_FILTER_LIST, ~ unique(.x[.x != ""]))
  ) %>%
  mutate(
    # Remove 'PASS' if there are other filters
    VARIANT_LEVEL_FILTER_LIST = map(
      VARIANT_LEVEL_FILTER_LIST, ~ if_else(
        length(.x) > 1 & "PASS" %in% .x,
        list(.x[.x != "PASS"]),
        list(.x)
      )[[1]]
    )
  ) %>%
  mutate(
    VARIANT_LEVEL_FILTER = map_chr(VARIANT_LEVEL_FILTER_LIST, ~ paste(.x, collapse = ";"))
  ) %>%
  select(
    -VARIANT_LEVEL_FILTER_STR_COMB,
    -VARIANT_LEVEL_FILTER_LIST
  ) %>%
  # Combine CHIP filters into a single semicolon-separated string
  mutate(
    CHIP_MUTATION_FILTER_STR_COMB = paste(
      CHIP_MUTATION_FILTER,
      PUTATIVE_CHIP_FILTER,
      sep = ";"
    )
  ) %>%
  mutate(
    CHIP_MUTATION_FILTER_LIST = str_split(CHIP_MUTATION_FILTER_STR_COMB, ";")
  ) %>%
  mutate(
    CHIP_MUTATION_FILTER_LIST = map(CHIP_MUTATION_FILTER_LIST, ~ unique(.x[.x != ""]))
  ) %>%
  mutate(
    CHIP_FILTER = map_chr(
      CHIP_MUTATION_FILTER_LIST, ~ case_when(
        length(.x) == 0 ~ "",
        length(.x) == 1 ~ .x[1],
        .default = paste(sort(.x), collapse = ";")
      )
    )
  ) %>%
  select(
    -CHIP_MUTATION_FILTER_STR_COMB,
    -CHIP_MUTATION_FILTER_LIST
  ) %>%
  # Create CHIP_INFO column
  mutate(
    CHIP_INFO = case_when(
      CHIP_MUTATION_FILTER == "" ~ paste0(
        "CHIP_Transcript=", Transcript, ";",
        "AAChange=", AAChange, ";",
        "GeneDetail=", GeneDetail, ";",
        "CHIP_Mutation_Class=", mutation_class, ";",
        "CHIP_Mutation_Definition=", mutation_definition, ";",
        "CHIP_Publication_Source=", gsub(";", "/", publication_source_concat, fixed = TRUE)
      ),
      .default = ""
    )
  ) %>%
  # Consolidate CHIP FILTER strings into one per variant
  group_by(
    CHROM,
    POS,
    REF,
    ALT
  ) %>%
  mutate(
    CHIP_VARIANT_LEVEL_FILTER = case_when(
      any(CHIP_FILTER == "") ~ "",
      .default = paste(CHIP_FILTER, collapse = ";")
    )
  ) %>%
  mutate(
    CHIP_VARIANT_LEVEL_FILTER_LIST = str_split(CHIP_VARIANT_LEVEL_FILTER, ";")
  ) %>%
  mutate(
    CHIP_VARIANT_LEVEL_FILTER_LIST = map(CHIP_VARIANT_LEVEL_FILTER_LIST, ~ unique(.x[.x != ""]))
  ) %>%
  mutate(
    CHIP_VARIANT_LEVEL_FILTER = map_chr(CHIP_VARIANT_LEVEL_FILTER_LIST, ~ paste(.x, collapse = ";"))
  ) %>%
  select(
    -CHIP_VARIANT_LEVEL_FILTER_LIST
  ) %>%
  # Consolidate CHIP INFO strings into one per variant
  mutate(
    CHIP_VARIANT_LEVEL_INFO = unique(CHIP_INFO[CHIP_INFO != ""]) %>% paste(collapse = ";")
  ) %>%
  ungroup() %>%
  # Create a combined filter column
  mutate(
    COMBINED_FILTER = case_when(
      VARIANT_LEVEL_FILTER == "PASS" & CHIP_FILTER == "" ~ "PASS",
      VARIANT_LEVEL_FILTER == "PASS" & CHIP_FILTER != "" ~ CHIP_FILTER,
      VARIANT_LEVEL_FILTER != "PASS" & CHIP_FILTER == "" ~ VARIANT_LEVEL_FILTER,
      .default = paste(VARIANT_LEVEL_FILTER, CHIP_FILTER, sep = ";") %>% str_replace_all("^;|;$", "")
    ),
    COMBINED_VARIANT_LEVEL_FILTER = case_when(
      VARIANT_LEVEL_FILTER == "PASS" & CHIP_VARIANT_LEVEL_FILTER == "" ~ "PASS",
      VARIANT_LEVEL_FILTER == "PASS" & CHIP_VARIANT_LEVEL_FILTER != "" ~ CHIP_VARIANT_LEVEL_FILTER,
      VARIANT_LEVEL_FILTER != "PASS" & CHIP_VARIANT_LEVEL_FILTER == "" ~ VARIANT_LEVEL_FILTER,
      .default = paste(VARIANT_LEVEL_FILTER, CHIP_VARIANT_LEVEL_FILTER, sep = ";") %>% str_replace_all("^;|;$", "")
    )
  )

# --- Add back in variants that failed to pass filters ---
df_filtered <- df_filtered %>%
  bind_rows(df_chip_tx_filtered) %>%
  # Handle NAs from bind
  mutate(
    SOMATICISM_FILTER = case_when(
      !is.na(SOMATICISM_FILTER) ~ SOMATICISM_FILTER,
      .default = ""
    ),
    EXONIC_SPLICING_FILTER = case_when(
      !is.na(EXONIC_SPLICING_FILTER) ~ EXONIC_SPLICING_FILTER,
      .default = ""
    ),
    CHIP_TRANSCRIPT_FILTER = case_when(
      !is.na(CHIP_TRANSCRIPT_FILTER) ~ CHIP_TRANSCRIPT_FILTER,
      .default = ""
    )
  ) %>%
  # Combine variant-level filters into a single semicolon-separated string
  group_by(
    CHROM,
    POS,
    REF,
    ALT
  ) %>%
  mutate(
    COMBINED_VARIANT_LEVEL_FILTER_STR_COMB = paste(
      FILTER,
      HARD_FILTER_AF,
      HARD_FILTER_VAF,
      HOMOPOLYMER_INDEL_FILTER,
      SOMATICISM_FILTER,
      EXONIC_SPLICING_FILTER,
      CHIP_TRANSCRIPT_FILTER,
      sep = ";"
    )
  ) %>%
  mutate(
    COMBINED_VARIANT_LEVEL_FILTER_STR_COMB = paste(COMBINED_VARIANT_LEVEL_FILTER_STR_COMB, collapse = ";")
  ) %>%
  ungroup() %>%
  mutate(
    COMBINED_VARIANT_LEVEL_FILTER_LIST = str_split(COMBINED_VARIANT_LEVEL_FILTER_STR_COMB, ";")
  ) %>%
  mutate(
    COMBINED_VARIANT_LEVEL_FILTER_LIST = map(COMBINED_VARIANT_LEVEL_FILTER_LIST, ~ unique(.x[.x != ""]))
  ) %>%
  mutate(
    # Remove 'PASS' (none of the variants have passed)
    COMBINED_VARIANT_LEVEL_FILTER_LIST = map(
      COMBINED_VARIANT_LEVEL_FILTER_LIST, ~ if_else(
        length(.x) > 1,
        list(.x[.x != "PASS"]),
        list(.x)
      )[[1]]
    )
  ) %>%
  mutate(
    COMBINED_VARIANT_LEVEL_FILTER = map_chr(COMBINED_VARIANT_LEVEL_FILTER_LIST, ~ paste(.x, collapse = ";"))
  ) %>%
  select(
    -COMBINED_VARIANT_LEVEL_FILTER_STR_COMB,
    -COMBINED_VARIANT_LEVEL_FILTER_LIST
  )

common_cols <- intersect(
  colnames(df),
  colnames(df_filtered)
)

df_final <- df %>%
  bind_rows(df_filtered[common_cols]) %>%
  select(
    CHROM,
    POS,
    ID,
    REF,
    ALT,
    QUAL,
    COMBINED_VARIANT_LEVEL_FILTER,
    CHIP_VARIANT_LEVEL_INFO,
    FORMAT,
    SAMPLE,
    Transcript,
    AAChange,
    GeneDetail,
    Func,
    ExonicFunc,
    gnomAD_AF,
    VAF,
    AD_REF,
    AD_ALT,
    F1R2_REF,
    F1R2_ALT,
    F2R1_REF,
    F2R1_ALT,
    SEQ_REF,
    SEQ_ALT,
    FILTER,
    HARD_FILTER_AF,
    HARD_FILTER_VAF,
    HOMOPOLYMER_FILTER,
    HOMOPOLYMER_INDEL_FILTER,
    SOMATICISM_FILTER,
    CHIP_MUTATION_FILTER,
    PUTATIVE_CHIP_FILTER,
    CHIP_INFO,
    all_of(colnames(chip_definitions)[colnames(chip_definitions) %in% colnames(df)])
  ) %>%
  rename(
    ORIGINAL_FILTER = FILTER,
    FILTER = COMBINED_VARIANT_LEVEL_FILTER,
    INFO = CHIP_VARIANT_LEVEL_INFO
  ) %>%
  # Remove NAs from FILTER and INFO columns
  mutate(
    FILTER = if_else(
      is.na(FILTER),
      ".",
      FILTER
    ),
    INFO = if_else(
      is.na(INFO),
      ".",
      INFO
    )
  ) %>%
  distinct() %>%
  # The COMBINED_VARIANT_LEVEL_FILTER and CHIP_VARIANT_LEVEL_INFO columns
  # *should* be the same for all rows with the same CHROM, POS, REF, and ALT
  # but in case they're not, we'll do one final merge
  group_by(
    CHROM,
    POS,
    REF,
    ALT
  ) %>%
  mutate(
    FILTER = paste(FILTER, collapse = ";"),
    INFO = unique(INFO[!(INFO %in% c("", "."))]) %>% paste(collapse = ";")
  ) %>%
  ungroup() %>%
  mutate(
    FILTER = str_split(FILTER, ";"),
    INFO = case_when(
      INFO == "" ~ ".",
      .default = INFO
    )
  ) %>%
  mutate(
    FILTER = map(FILTER, ~ unique(.x[.x != ""])),
  ) %>%
  mutate(
    # Remove 'PASS' if there are other filters
    FILTER = map(
      FILTER, ~ if_else(
        length(.x) > 1 & "PASS" %in% .x,
        list(.x[.x != "PASS"]),
        list(.x)
      )[[1]]
    )
  ) %>%
  mutate(
    FILTER = map_chr(FILTER, ~ paste(.x, collapse = ";")),
  ) %>%
  distinct() %>%
  # Merge INFO fields
  # For example, if a variant has multiple CHIP transcripts, merge the CHIP_Transcript field:
  # CHIP_Transcript=NM_0000001|NM_0000002
  mutate(
    INFO_split = str_split(INFO, ";")
  ) %>%
  mutate(
    INFO_CHIP_Transcript = map_chr(
      INFO_split, ~ (
        grep("^CHIP_Transcript=", .x, value = TRUE) %>%
          gsub("^CHIP_Transcript=", "", .) %>%
          paste(collapse = "|")
      )
    ),
    INFO_AAChange = map_chr(
      INFO_split, ~ (
        grep("^AAChange=", .x, value = TRUE) %>%
          gsub("^AAChange=", "", .) %>%
          paste(collapse = "|")
      )
    ),
    INFO_GeneDetail = map_chr(
      INFO_split, ~ (
        grep("^GeneDetail=", .x, value = TRUE) %>%
          gsub("^GeneDetail=", "", .) %>%
          paste(collapse = "|")
      )
    ),
    INFO_CHIP_Mutation_Class = map_chr(
      INFO_split, ~ (
        grep("^CHIP_Mutation_Class=", .x, value = TRUE) %>%
          gsub("^CHIP_Mutation_Class=", "", .) %>%
          paste(collapse = "|")
      )
    ),
    INFO_CHIP_Mutation_Definition = map_chr(
      INFO_split, ~ (
        grep("^CHIP_Mutation_Definition=", .x, value = TRUE) %>%
          gsub("^CHIP_Mutation_Definition=", "", .) %>%
          paste(collapse = "|")
      )
    ),
    INFO_CHIP_Publication_Source = map_chr(
      INFO_split, ~ (
        grep("^CHIP_Publication_Source=", .x, value = TRUE) %>%
          gsub("^CHIP_Publication_Source=", "", .) %>%
          paste(collapse = "|")
      )
    )
  ) %>%
  mutate(
    INFO = case_when(
      INFO_CHIP_Transcript == "" ~ "CHIP_Transcript=.;AAChange=.;GeneDetail=.;CHIP_Mutation_Class=.;CHIP_Mutation_Definition=.;CHIP_Publication_Source=.",
      .default = paste0(
        "CHIP_Transcript=", INFO_CHIP_Transcript, ";",
        "AAChange=", INFO_AAChange, ";",
        "GeneDetail=", INFO_GeneDetail, ";",
        "CHIP_Mutation_Class=", INFO_CHIP_Mutation_Class, ";",
        "CHIP_Mutation_Definition=", INFO_CHIP_Mutation_Definition, ";",
        "CHIP_Publication_Source=", INFO_CHIP_Publication_Source
      )
    )
  ) %>%
  select(
    -INFO_split,
    -INFO_CHIP_Transcript,
    -INFO_AAChange,
    -INFO_GeneDetail,
    -INFO_CHIP_Mutation_Class,
    -INFO_CHIP_Mutation_Definition,
    -INFO_CHIP_Publication_Source
  ) %>%
  distinct()

# --- Construct output VCF for annotation ---
df_final_vcf <- df_final %>%
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
  distinct() %>%
  group_by(
    CHROM,
    POS,
    REF
  ) %>%
  mutate(
    INFO_CHIP_Multiallelic_Filters = case_when(
      length(unique(FILTER)) > 1 ~ paste0("CHIP_Multiallelic_Filters=", gsub(";", "|", FILTER, fixed = TRUE)),
      .default = ""
    )
  ) %>%
  ungroup() %>%
  mutate(
    INFO = case_when(
      INFO_CHIP_Multiallelic_Filters != "" ~ paste0(INFO, ";", INFO_CHIP_Multiallelic_Filters),
      .default = INFO
    )
  ) %>%
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

# --- Write output VCF ---
# Read in header from input VCF
vcf_header <- read_lines(args$vcf_header) %>%
  keep(~ str_detect(.x, "^#"))

# Define new FILTER fields
filter_fields <- list(
  gnomad_af_na = "Variant is not present in gnomAD, and NAs were not allowed to be interpreted as 0",
  chip_gnomad_af_filter_fail = "Variant has a gnomAD AF >= 0.01",
  chip_vaf_filter_fail = "Variant has a VAF < 0.02",
  homopolymer_variant = "Variant is within a homopolymer region",
  chip_homopolymer_indel_filter_fail = "Variant is an INDEL within a homopolymer region and has AD_ALT < 10 or VAF < 0.1",
  not_exonic_or_splicing_variant = "Variant is not exonic or splicing",
  chip_somaticism_filter_fail = "Variant is within a transcript specified in the somaticism_transcripts file and did not pass a binomial test (n = AD_ALT, k = DP, p = 0.5, alpha = 0.001)",
  exonic_or_splicing_variant_not_in_chip_transcript = "Variant is exonic or splicing but is not within a CHIP transcript",
  mutation_in_c_terminal = "Variant is within a CHIP gene but is in the C-terminal domain of the protein",
  chip_mutation_match_filter_fail = "Variant does not match a CHIP mutation definition",
  chip_putative_filter_fail = "Variant is putative and has AD_ALT < 5, VAF > 0.2, F1R2_REF < 2, F1R2_ALT < 2, F2R1_REF < 2, or F2R1_ALT < 2"
)

create_filter_header <- function(filter_name, filter_description) {
  paste0(
    "##FILTER=<ID=", filter_name, ",Description=\"", filter_description, "\">"
  )
}

filter_fields_header <- map_chr(
  names(filter_fields),
  ~ create_filter_header(.x, filter_fields[[.x]])
)

info_fields <- list(
  CHIP_Transcript = "Transcript ID of the variant",
  AAChange = "Amino acid change of the variant",
  GeneDetail = "Gene detail of the variant",
  CHIP_Mutation_Class = "Mutation class of the variant",
  CHIP_Mutation_Definition = "Mutation definition of the variant",
  CHIP_Publication_Source = "Publication source of the mutation definition",
  CHIP_Multiallelic_Filters = "Allele-specific filters for multiallelic variants"
)

create_info_header <- function(info_name, info_description) {
  paste0(
    "##INFO=<ID=", info_name, ",Number=A,Type=String,Description=\"", info_description, "\">"
  )
}

info_fields_header = map_chr(
  names(info_fields),
  ~ create_info_header(.x, info_fields[[.x]])
)

# Insert new FILTER and INFO fields into header
vcf_header_filter_end <- grep("^##FILTER", vcf_header)
if (length(vcf_header_filter_end) == 0) {
  vcf_header_filter_end <- grep("^#CHROM", vcf_header) %>% max() - 1
} else {
  vcf_header_filter_end <- vcf_header_filter_end %>% max()
}
vcf_header_info_end <- grep("^##INFO", vcf_header)
if (length(vcf_header_info_end) == 0) {
  vcf_header_info_end <- grep("^#CHROM", vcf_header) %>% max() - 1
} else {
  vcf_header_info_end <- vcf_header_info_end %>% max()
}
filter_first <- vcf_header_filter_end < vcf_header_info_end
vcf_header_length <- length(vcf_header)

if (filter_first) {
  vcf_header_1 <- vcf_header[1:vcf_header_filter_end]
  vcf_header_2 <- vcf_header[(vcf_header_filter_end + 1):vcf_header_info_end]
  new_vcf_header <- c(vcf_header_1, filter_fields_header, vcf_header_2, info_fields_header)
  if (vcf_header_length > vcf_header_info_end) {
    new_vcf_header <- c(new_vcf_header, vcf_header[(vcf_header_info_end + 1):vcf_header_length])
  }
} else {
  vcf_header_1 <- vcf_header[1:vcf_header_info_end]
  vcf_header_2 <- vcf_header[(vcf_header_info_end + 1):vcf_header_filter_end]
  new_vcf_header <- c(vcf_header_1, info_fields_header, vcf_header_2, filter_fields_header)
  if (vcf_header_length > vcf_header_filter_end) {
    new_vcf_header <- c(new_vcf_header, vcf_header[(vcf_header_filter_end + 1):vcf_header_length])
  }
}

# Create a copy of the header with the FORMAT and SAMPLE columns removed from the last line
new_vcf_header_no_format <- new_vcf_header
new_vcf_header_no_format[length(new_vcf_header_no_format)] <- gsub("\\tFORMAT.*$", "", new_vcf_header_no_format[length(new_vcf_header_no_format)])

# Write header to output VCF
output_vcf <- paste0(args$output_prefix, ".vcf")
write_lines(new_vcf_header_no_format, output_vcf)

# Write variants to output VCF
write_tsv(df_final_vcf, output_vcf, append = TRUE, col_names = FALSE)

# --- Write output CSVs ---
output_csv <- paste0(args$output_prefix, ".csv")
# Replace "." and NAs in all columns with empty strings
df_final_csv <- df_final %>%
  # First, convert integer, double, and logical columns to characters
  # Ensure that numbers are not converted to scientific notation
  mutate(
    across(where(is.integer), ~ format(.x, scientific = FALSE, trim = TRUE, justify = "none")),
    across(where(is.double), ~ format(.x, scientific = FALSE, trim = TRUE, justify = "none")),
    across(where(is.logical), ~ if_else(.x, "TRUE", "FALSE"))
  ) %>%
  mutate(
    across(everything(), ~ if_else(is.na(.), "", .)),
    across(everything(), ~ if_else(. == ".", "", .))
  )
write_csv(df_final_csv, output_csv)

# --- Write RData ---
save.image(paste0(args$output_prefix, ".RData"))