apply_hard_filters <- function(df, tumor_sample_name) {
  # Filter by AD, DP, AF, F1R2/F2R1 VCF fields and gnomAD frequency
  get_format_field <- function(x, format_field) { grep(paste("^", x, "$", sep = ""), format_field, perl = TRUE) }
  # vars_stats <- apply(df[c("FORMAT", tumor_sample_name, "gnomAD_AF")], 1, function(x) {
  vars_stats <- apply(df[c("FORMAT", tumor_sample_name, "gnomAD_AF")], 1, function(x) {
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
    var_n_vaf <- (var_n_ad / var_n_dp)[2:length(var_n_ad)]
    var_n_vaf_gte_2pc <- var_n_vaf >= 0.02
    var_n_vaf_lt_35pc <- var_n_vaf < 0.35
    var_n_vaf_m2_i <- get_format_field("AF", format_field)
    var_n_vaf_m2 <- as.numeric(strsplit(sample_field[var_n_vaf_m2_i], ",")[[1]])
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
    HARD_FILTER = "FAIL"
    if(filter_pass) {
      HARD_FILTER = "PASS"
    } else if(filter_manual_review) {
      HARD_FILTER = "MANUAL_REVIEW"
    }
    return(data.frame(
      AD = paste(as.character(var_n_ad), collapse = ","),
      DP = as.character(var_n_dp),
      VAF = paste(as.character(var_n_vaf), collapse = ","),
      VAF_M2 = paste(as.character(var_n_vaf_m2), collapse = ","),
      F1R2 = paste(as.character(var_n_for), collapse = ","),
      F2R1 = paste(as.character(var_n_rev), collapse = ","),
      HARD_FILTER = HARD_FILTER
    ))
  })
  vars_stats <- do.call(rbind, vars_stats)
  df <- cbind(df, vars_stats)
  df$COMBINED_FILTER = apply(df[c("FILTER", "HARD_FILTER")], 1, function(x) {
    if(x[[1]] == "PASS") {
      return(x[[2]])
    } else {
      return("FAIL")
    }
  })
  return(df)
}