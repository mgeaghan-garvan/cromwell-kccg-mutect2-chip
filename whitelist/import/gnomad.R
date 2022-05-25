rename_gnomad_col <- function(df, gnomad_source) {
    # Ensure the following columns are present
    req_cols <- c("Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene", "GeneDetail.refGene", "ExonicFunc.refGene", "AAChange.refGene")
    gnomad_g_cols <- NA
    gnomad_e_cols <- NA
    gnomad_cols <- NA
    new_gnomad_g_cols <- NA
    new_gnomad_e_cols <- NA
    new_gnomad_cols <- NA
    # When using both genome and exome gnomAD annotations, Annovar (unhelpfully) gives them the same column headers, with the second annotation suffixed with ".1".
    # So, we need to supply the script with the order in which they were annotated ("genome,exome", or "exome,genome") so it can determine which is which.
    if (gnomad_source == "genome,exome") {
        gnomad_g_cols <- c("AF", "AF_popmax", "AF_male", "AF_female", "AF_raw", "AF_afr", "AF_sas", "AF_amr", "AF_eas", "AF_nfe", "AF_fin", "AF_asj", "AF_oth", "non_topmed_AF_popmax", "non_neuro_AF_popmax", "non_cancer_AF_popmax", "controls_AF_popmax")
        gnomad_e_cols <- c("AF.1", "AF_popmax.1", "AF_male.1", "AF_female.1", "AF_raw.1", "AF_afr.1", "AF_sas.1", "AF_amr.1", "AF_eas.1", "AF_nfe.1", "AF_fin.1", "AF_asj.1", "AF_oth.1", "non_topmed_AF_popmax.1", "non_neuro_AF_popmax.1", "non_cancer_AF_popmax.1", "controls_AF_popmax.1")
        req_cols <- c(req_cols, gnomad_g_cols, gnomad_e_cols)
        new_gnomad_g_cols <- paste("gnomAD_genome_", gnomad_g_cols, sep = "")
        new_gnomad_e_cols <- gsub("\\.\\d+$", "", gnomad_e_cols, perl = TRUE)
        new_gnomad_e_cols <- paste("gnomAD_exome_", new_gnomad_e_cols, sep = "")
    } else if (gnomad_source == "exome,genome") {
        gnomad_e_cols <- c("AF", "AF_popmax", "AF_male", "AF_female", "AF_raw", "AF_afr", "AF_sas", "AF_amr", "AF_eas", "AF_nfe", "AF_fin", "AF_asj", "AF_oth", "non_topmed_AF_popmax", "non_neuro_AF_popmax", "non_cancer_AF_popmax", "controls_AF_popmax")
        gnomad_g_cols <- c("AF.1", "AF_popmax.1", "AF_male.1", "AF_female.1", "AF_raw.1", "AF_afr.1", "AF_sas.1", "AF_amr.1", "AF_eas.1", "AF_nfe.1", "AF_fin.1", "AF_asj.1", "AF_oth.1", "non_topmed_AF_popmax.1", "non_neuro_AF_popmax.1", "non_cancer_AF_popmax.1", "controls_AF_popmax.1")
        req_cols <- c(req_cols, gnomad_g_cols, gnomad_e_cols)
        new_gnomad_e_cols <- paste("gnomAD_exome_", gnomad_e_cols, sep = "")
        new_gnomad_g_cols <- gsub("\\.\\d+$", "", gnomad_g_cols, perl = TRUE)
        new_gnomad_g_cols <- paste("gnomAD_genome_", new_gnomad_g_cols, sep = "")
    } else if (gnomad_source == "exome" || gnomad_source == "genome") {
        gnomad_cols <- c("AF", "AF_popmax", "AF_male", "AF_female", "AF_raw", "AF_afr", "AF_sas", "AF_amr", "AF_eas", "AF_nfe", "AF_fin", "AF_asj", "AF_oth", "non_topmed_AF_popmax", "non_neuro_AF_popmax", "non_cancer_AF_popmax", "controls_AF_popmax")
        req_cols <- c(req_cols, gnomad_cols)
        new_gnomad_cols <- paste("gnomAD_", gnomad_source, "_", gnomad_cols, sep = "")
    } else {
        stop("Invalid gnomAD source.")
    }
    stopifnot(all(req_cols %in% colnames(df)))
    if (gnomad_source %in% c("genome,exome", "exome,genome")) {
        colnames(df)[colnames(df) %in% gnomad_g_cols] <- new_gnomad_g_cols
        colnames(df)[colnames(df) %in% gnomad_e_cols] <- new_gnomad_e_cols
    } else {
        colnames(df)[colnames(df) %in% gnomad_cols] <- new_gnomad_cols
    }
    return(df)
}

get_gnomad_af <- function(df, gnomad_source, gnomad_pop, treat_missing_as_rare) {
    # Get gnomAD allele frequencies
    if (gnomad_source %in% c("genome,exome", "exome,genome")) {
        gnomad_pop_columns <- paste("gnomAD", c("genome", "exome"), gnomad_pop, sep = "_")

        df$gnomAD_AF <- apply(df[gnomad_pop_columns], 1, function(x) {
            af_g <- x[[1]]
            af_e <- x[[2]]
            missing_vals <- c(".", "", NA)
            if(!(af_e %in% missing_vals)) {
                return(as.numeric(af_e))
            } else if(!(af_g %in% missing_vals)) {
                return(as.numeric(af_g))
            } else if(treat_missing_as_rare) {
                return(0)
            } else {
                return(NA)
            }
        })
    } else {
        gnomad_pop_columns <- paste("gnomAD", gnomad_source, gnomad_pop, sep = "_")

        df$gnomAD_AF <- apply(df[gnomad_pop_columns], 1, function(x) {
            af <- x[[1]]
            missing_vals <- c(".", "", NA)
            if(!(af %in% missing_vals)) {
                return(as.numeric(af))
            } else if(treat_missing_as_rare) {
                return(0)
            } else {
                return(NA)
            }
        })
    }
    return(df)
}