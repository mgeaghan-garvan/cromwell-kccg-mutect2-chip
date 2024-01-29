#!/bin/bash

INPUT_VCF=${1}
CHIP_MUTATIONS_CSV=${2}
REF_FASTA=${3}
SOMATICISM_TRANSCRIPTS=${4}

set -euo pipefail

# === STEP 1: Create BED file of CHIP gene regions ===

CHIP_MUTATIONS_TMP_BED="chip_mutations.genomic_regions.bed"
CHIP_MUTATIONS_BED="chip_mutations.genomic_regions.merged.bed"

R --vanilla <<EOF
library(tidyverse)
df <- read_csv("${CHIP_MUTATIONS_CSV}")
df <- df %>% select(c(chr, gene_genomic_start, gene_genomic_end)) %>% mutate(gene_genomic_start = gene_genomic_start - 1) %>% distinct %>% arrange(chr, gene_genomic_start)
write_tsv(df, "${CHIP_MUTATIONS_TMP_BED}", col_names = FALSE)
EOF

# Merge overlapping regions
bedtools merge -i ${CHIP_MUTATIONS_TMP_BED} > ${CHIP_MUTATIONS_BED}

# === STEP 2: Filter input VCF for CHIP gene regions ===
# INPUT_VCF_BN="$(echo ${INPUT_VCF} | sed -E -e 's/\.vcf(\.gz)?$//g')"
INPUT_VCF_BN="$(basename $(basename ${INPUT_VCF} .gz) .vcf)"
CHIP_GENES_VCF="${INPUT_VCF_BN}.chip_genes.vcf"
NON_CHIP_GENES_SITES_ONLY_VCF="${INPUT_VCF_BN}.non_chip_genes.so.vcf"
bedtools intersect -a ${INPUT_VCF} -b ${CHIP_MUTATIONS_BED} -wa -header > ${CHIP_GENES_VCF}

# === STEP 3: Create a non-CHIP gene regions VCF ===
bedtools intersect -a ${INPUT_VCF} -b ${CHIP_MUTATIONS_BED} -wa -header -v | cut -f 1-8 > ${NON_CHIP_GENES_SITES_ONLY_VCF}

# Update FILTER column and strip the INFO column
NON_CHIP_GENES_SITES_ONLY_FILTER_VCF="${INPUT_VCF_BN}.non_chip_genes.so.filter.vcf"
awk -v FS="\t" -v OFS="\t" '
	$0 ~ /^#/ { print $0 }
	$0 !~ /^#/ {
	print $1, $2, $3, $4, $5, $6, "non_chip_gene", "."
	}
' ${NON_CHIP_GENES_SITES_ONLY_VCF} > ${NON_CHIP_GENES_SITES_ONLY_FILTER_VCF}

# Update header with new FILTER information
NEW_FILTER_INFO="##FILTER=<ID=non_chip_gene,Description=\"Variant is not in a CHIP gene region\">"
# Find line number of last FILTER line
LAST_FILTER_LINE_NUM=$(grep -n "^##FILTER" ${NON_CHIP_GENES_SITES_ONLY_FILTER_VCF} | tail -n 1 | cut -d ":" -f 1)
# Insert new FILTER line after last FILTER line
sed -i "${LAST_FILTER_LINE_NUM}a ${NEW_FILTER_INFO}" ${NON_CHIP_GENES_SITES_ONLY_FILTER_VCF}

# BGZIP and tabix the non-CHIP gene annotations
bgzip -c ${NON_CHIP_GENES_SITES_ONLY_FILTER_VCF} > ${NON_CHIP_GENES_SITES_ONLY_FILTER_VCF}.gz
tabix -f -s 1 -b 2 -e 2 ${NON_CHIP_GENES_SITES_ONLY_FILTER_VCF}.gz
NON_CHIP_GENES_SITES_ONLY_FILTER_VCF="${NON_CHIP_GENES_SITES_ONLY_FILTER_VCF}.gz"

# === STEP 4: Create normalized VCF with multi-allelic variants split ===
CHIP_GENES_VCF_BN="$(basename ${CHIP_GENES_VCF} .vcf)"
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
    -e 'FORMAT/F1R2[0:*]<1' | \
bcftools filter \
    -s chip_f2r1_filter_fail \
    -m + \
    -e 'FORMAT/F2R1[0:*]<1' > ${FILTERED_CHIP_GENES_NORM_NO_INFO_VCF}

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

# === STEP 10: Get the header from the VCF ===
FILTERED_CHIP_GENES_NORM_NO_INFO_VCF_HEADER="${FILTERED_CHIP_GENES_NORM_NO_INFO_VCF}.header"
grep "^#" ${FILTERED_CHIP_GENES_NORM_NO_INFO_VCF} > ${FILTERED_CHIP_GENES_NORM_NO_INFO_VCF_HEADER}

# === STEP 10: Run R script ===
CHIP_ANNOTATION_PREFIX="${CHIP_GENES_VCF_BN}.norm.no_info.hard_filter.chip_annotations"
Rscript --vanilla chip_annotation.R \
    --sample ${INPUT_VCF_BN} \
    --vcf_header ${FILTERED_CHIP_GENES_NORM_NO_INFO_VCF_HEADER} \
    --chip_definitions ${CHIP_MUTATIONS_CSV} \
    --seq ${SEQ_TSV} \
    --annovar ${ANNOVAR_PREFIX}.hg38_multianno.txt \
    --annovar_function ${ANNOVAR_PREFIX}.refGene.variant_function \
    --annovar_exonic_function ${ANNOVAR_PREFIX}.refGene.exonic_variant_function \
    --somaticism_transcripts ${SOMATICISM_TRANSCRIPTS} \
    --output ${CHIP_ANNOTATION_PREFIX}

# BGZIP and tabix the CHIP annotations
bcftools sort ${CHIP_ANNOTATION_PREFIX}.vcf > ${CHIP_ANNOTATION_PREFIX}.sorted.vcf
bgzip -c ${CHIP_ANNOTATION_PREFIX}.sorted.vcf > ${CHIP_ANNOTATION_PREFIX}.sorted.vcf.gz
tabix -f -s 1 -b 2 -e 2 ${CHIP_ANNOTATION_PREFIX}.sorted.vcf.gz

# === STEP 11: Annotate the original VCF with both CHIP annotations and non-CHIP gene annotations ===
TMP_VCF="${INPUT_VCF_BN}.chip.tmp.vcf"
FINAL_VCF="${INPUT_VCF_BN}.chip.vcf"
bcftools annotate \
    -a ${NON_CHIP_GENES_SITES_ONLY_FILTER_VCF} \
    -c "=FILTER" \
    -o ${TMP_VCF} \
    ${INPUT_VCF}
bgzip -c ${TMP_VCF} > ${TMP_VCF}.gz
tabix -f -s 1 -b 2 -e 2 ${TMP_VCF}.gz
bcftools annotate \
    -a ${CHIP_ANNOTATION_PREFIX}.sorted.vcf.gz \
    -c "=FILTER,=INFO" \
    -o ${FINAL_VCF} \
    ${TMP_VCF}.gz
bgzip -c ${FINAL_VCF} > ${FINAL_VCF}.gz
tabix -f -s 1 -b 2 -e 2 ${FINAL_VCF}.gz

rm ${TMP_VCF} ${TMP_VCF}.gz ${TMP_VCF}.gz.tbi

echo DONE