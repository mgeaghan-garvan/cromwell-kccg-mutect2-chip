apply_chip_filters <- function(df, match_ns_func, match_fs_func, match_sg_func, match_sp_func, parse_cond_func, parse_aa_change_func) {
    # Matching non-synonymous mutations (SNV, non-frameshift INDELs)
    wl_ns <- apply(vars_g_chip_func_filtered[, c("ExonicFunc.refGene", "Nonsynonymous", "Putative_Nonsynonymous", "AAChange.protChange", "AAChange.exon")], 1, match_ns, parse_cond_func = parse_cond, parse_aa_change_func = parse_aa_change)

    # Matching frameshift mutations
    wl_fs <- apply(vars_g_chip_func_filtered[, c("ExonicFunc.refGene", "Frameshift", "Putative_Frameshift", "AAChange.protChange", "AAChange.exon")], 1, match_fs, parse_cond_func = parse_cond, parse_aa_change_func = parse_aa_change)

    # Matching stop-gain mutations
    wl_sg <- apply(vars_g_chip_func_filtered[, c("ExonicFunc.refGene", "Stop_Gain", "Putative_Stop_Gain", "AAChange.protChange", "AAChange.exon")], 1, match_sg, parse_cond_func = parse_cond, parse_aa_change_func = parse_aa_change)

    # Matching splicing mutations
    wl_sp <- apply(vars_g_chip_func_filtered[, c("Func.refGene", "Splicing", "Putative_Splicing", "GeneDetail.transcript", "GeneDetail.exon")], 1, match_sp, parse_cond_func = parse_cond)

    # Add whitelist/manual review columns to data frame
    wl_ns_df <- data.frame(t(wl_ns))
    colnames(wl_ns_df) <- paste(c("Whitelist", "Manual_Review", "Putative_Whitelist", "Putative_Manual_Review"), "Nonsynonymous", sep = "_")

    wl_fs_df <- data.frame(t(wl_fs))
    colnames(wl_fs_df) <- paste(c("Whitelist", "Manual_Review", "Putative_Whitelist", "Putative_Manual_Review"), "Frameshift", sep = "_")

    wl_sg_df <- data.frame(t(wl_sg))
    colnames(wl_sg_df) <- paste(c("Whitelist", "Manual_Review", "Putative_Whitelist", "Putative_Manual_Review"), "Stop_Gain", sep = "_")

    wl_sp_df <- data.frame(t(wl_sp))
    colnames(wl_sp_df) <- paste(c("Whitelist", "Manual_Review", "Putative_Whitelist", "Putative_Manual_Review"), "Splicing", sep = "_")

    # Overall whitelist/manual review
    wl <- wl_ns_df[, 1] | wl_fs_df[, 1] | wl_sg_df[, 1] | wl_sp_df[, 1]
    mr <- (!wl) & (wl_ns_df[, 2] | wl_fs_df[, 2] | wl_sg_df[, 2] | wl_sp_df[, 2])
    pwl <- wl_ns_df[, 3] | wl_fs_df[, 3] | wl_sg_df[, 3] | wl_sp_df[, 3]
    pmr <- (!pwl) & (wl_ns_df[, 4] | wl_fs_df[, 4] | wl_sg_df[, 4] | wl_sp_df[, 4])
    overall_wl <- data.frame(Whitelist = wl, Manual_Review = mr, Putative_Whitelist = pwl, Putative_Manual_Review = pmr)
    
    whitelist <- list()
    whitelist$wl_ns_df <- wl_ns_df
    whitelist$wl_fs_df <- wl_fs_df
    whitelist$wl_sg_df <- wl_sg_df
    whitelist$wl_sp_df <- wl_sp_df
    whitelist$overall_wl <- overall_wl
    return(whitelist)
}