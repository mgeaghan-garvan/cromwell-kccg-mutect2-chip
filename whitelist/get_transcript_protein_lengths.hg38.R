library(biomaRt)

ensembl <- useMart("ENSEMBL_MART_ENSEMBL", "hsapiens_gene_ensembl")

refseq <- getBM(attributes = c("refseq_mrna", "ensembl_transcript_id", "hgnc_symbol", "ensembl_gene_id"), mart = ensembl)

cds <- getBM(attributes = c("ensembl_transcript_id", "cds_length"), mart = ensembl)
refseq_cds <- merge(refseq, cds)
refseq_cds <- unique(refseq_cds[c("refseq_mrna", "hgnc_symbol", "cds_length")])
refseq_cds$prot_length <- as.integer(refseq_cds$cds_length / 3)
refseq_cds <- refseq_cds[!is.na(refseq_cds$prot_length),]
refseq_cds_valid <- refseq_cds[refseq_cds$refseq_mrna != "",]

write.csv(refseq_cds, "transcript_protein_lengths.all.csv", row.names = FALSE)
write.csv(refseq_cds_valid, "transcript_protein_lengths.csv", row.names = FALSE)
