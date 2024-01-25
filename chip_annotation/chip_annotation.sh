#!/bin/bash

INPUT_VCF=${1}
REF_FASTA=${2}

set -euo pipefail

# === STEP 1: Create BED file of CHIP gene regions ===

CHIP_MUTATIONS_CSV="chip_mutations/chip_mutations.csv"
CHIP_MUTATIONS_BED="chip_mutations.genomic_regions.merged.bed"

R --vanilla <<EOF
library(tidyverse)
df <- read_csv("chip_mutations/chip_mutations.csv")
df <- df %>% select(c(chr, gene_genomic_start, gene_genomic_end)) %>% mutate(gene_genomic_start = gene_genomic_start - 1) %>% distinct %>% arrange(chr, gene_genomic_start)
write_tsv(df, "chip_mutations.genomic_regions.bed", col_names = FALSE)
EOF

# Merge overlapping regions
bedtools merge -i chip_mutations.genomic_regions.bed > ${CHIP_MUTATIONS_BED}

# === STEP 2: Filter input VCF for CHIP gene regions ===
INPUT_VCF_BN="$(echo ${INPUT_VCF} | sed -E -e 's/\.vcf(\.gz)?$//g')"
CHIP_GENES_VCF="${INPUT_VCF_BN}.chip_genes.vcf"
NON_CHIP_GENES_SITES_ONLY_VCF="${INPUT_VCF_BN}.non_chip_genes.so.vcf"
bedtools intersect -a ${INPUT_VCF} -b ${CHIP_MUTATIONS_BED} -wa -header > ${CHIP_GENES_VCF}

# === STEP 3: Create a non-CHIP gene regions VCF ===
bedtools intersect -a ${INPUT_VCF} -b ${CHIP_MUTATIONS_BED} -wa -header -v | cut -f 1-8 > ${NON_CHIP_GENES_SITES_ONLY_VCF}

# Append FILTER and INFO columns
NON_CHIP_GENES_SITES_ONLY_FILTER_VCF="${INPUT_VCF_BN}.non_chip_genes.so.filter.vcf"
awk -v FS="\t" -v OFS="\t" '
	$0 ~ /^#/ { print $0 }
	$0 !~ /^#/ {
	if (($7 == "PASS") || ($7 == "") || ($7 == ".")) { p = "non_chip_gene" } else { p = $7";non_chip_gene" };
	if (($8 == "") || ($8 == ".") ) { i = "CHIP_FILTER=skipped_non_chip_gene" } else { i = $8";CHIP_FILTER=skipped_non_chip_gene" };
	print $1, $2, $3, $4, $5, $6, p, i
	}
' ${NON_CHIP_GENES_SITES_ONLY_VCF} > ${NON_CHIP_GENES_SITES_ONLY_FILTER_VCF}

# === STEP 4: Create normalized VCF with multi-allelic variants split ===
CHIP_GENES_VCF_BN="$(echo ${CHIP_GENES_VCF} | sed -E -e 's/\.vcf(\.gz)?$//g')"
CHIP_GENES_NORM_VCF="${CHIP_GENES_VCF_BN}.norm.vcf"
bcftools norm -m -any -o ${CHIP_GENES_NORM_VCF} ${CHIP_GENES_VCF}

# === STEP 5: Create BED file with 10bp flanking regions ===

CHIP_GENES_NORM_BED="${INPUT_VCF_BN}.norm.bed"
awk -v OFS="\t" '$0 !~ /^#/ {print $1, $2 - 11, $2 + length($4) + 9, $1":"$2":"$4":"$5}' ${CHIP_GENES_NORM_VCF} > ${CHIP_GENES_NORM_BED}

# === STEP 6: Create a sequence context TSV ===
SEQ_TSV="${INPUT_VCF_BN}.seq.tsv"
bedtools getfasta -fi ${REF_FASTA} -bed ${CHIP_GENES_NORM_BED} -tab -nameOnly > ${SEQ_TSV}

# === STEP 7: Strip INFO field from CHIP gene regions VCF ===
CHIP_GENES_NORM_NO_INFO_VCF="${CHIP_GENES_VCF_BN}.norm.no_info.vcf"
bcftools annotate -x INFO -o ${CHIP_GENES_NORM_NO_INFO_VCF} ${CHIP_GENES_NORM_VCF}

# === STEP 8: Apply hard filters to CHIP genes VCF ===
FILTERED_CHIP_GENES_NORM_NO_INFO_VCF="${CHIP_GENES_VCF_BN}.norm.no_info.hard_filter.vcf"
cat ${CHIP_GENES_NORM_NO_INFO_VCF} | \
bcftools filter \
    -s chip_ad_filter_fail \
    -m + \
    -i 'FORMAT/AD[0:1]>=3' | \
bcftools filter \
    -s chip_dp_filter_fail \
    -m + \
    -i 'FORMAT/DP[0]>=20' | \
bcftools filter \
    -s chip_f1r2_filter_fail \
    -m + \
    -i 'FORMAT/F1R2[0:0]>=1' | \
bcftools filter \
    -s chip_f1r2_filter_fail \
    -m + \
    -i 'FORMAT/F1R2[0:1]>=1' | \
bcftools filter \
    -s chip_f2r1_filter_fail \
    -m + \
    -i 'FORMAT/F2R1[0:0]>=1' | \
bcftools filter \
    -s chip_f2r1_filter_fail \
    -m + \
    -i 'FORMAT/F2R1[0:1]>=1' > ${FILTERED_CHIP_GENES_NORM_NO_INFO_VCF}

# === STEP 9: Run Annovar on filtered VCF ===
ANNOVAR_PREFIX="${CHIP_GENES_VCF_BN}.norm.no_info.hard_filter.annot"
table_annovar.pl \
    ${FILTERED_CHIP_GENES_NORM_NO_INFO_VCF} \
    ~/apps/annovar/humandb \
    -buildver hg38 \
    -out ${ANNOVAR_PREFIX} \
    -protocol refGene,ensGene,gnomad211_genome,gnomad211_exome \
    -operation g,g,f,f \
    -nastring . \
    -vcfinput \
    -polish \
    -argument "-exonicsplicing -transcript_function -separate,-exonicsplicing -transcript_function -separate,,"

# === STEP 10: Run R script ===

echo DONE