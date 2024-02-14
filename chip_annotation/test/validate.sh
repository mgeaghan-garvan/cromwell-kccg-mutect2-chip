#!/bin/bash
# Validate test outputs

OUTDIR="output"
POS_VCF="${OUTDIR}/positive_sites.chip.vcf.gz"
NEG_VCF="${OUTDIR}/negative_sites.chip.vcf.gz"
ALL_VCF="${OUTDIR}/all_sites.chip.vcf.gz"

# === Positive tests ===
# Test that all FILTER values are PASS
if [ "$(zgrep -v '^#' ${POS_VCF} | awk '$7 != "PASS"' | wc -l)" -ne 0 ]; then
    echo "FAIL: Positive sites VCF contains non-PASS FILTER values"
    exit 1
fi

# Test that no INFO fields contain an empty CHIP mutation match
# i.e. "CHIP_Transcript=.;" or "CHIP_Transcript=.,.;"
if [ "$(zgrep -v '^#' ${POS_VCF} | awk '$8 ~ /CHIP_Transcript=\.(,\.)*;/' | wc -l)" -ne 0 ]; then
    echo "FAIL: Positive sites VCF contains empty CHIP mutation matches in INFO field"
    exit 1
fi

# === Negative tests ===
# Test that no FILTER values are PASS
if [ "$(zgrep -v '^#' ${NEG_VCF} | awk '$7 == "PASS"' | wc -l)" -ne 0 ]; then
    echo "FAIL: Negative sites VCF contains PASS FILTER values"
    exit 1
fi

# Test that all sites that fail due to a CHIP mutation mismatch
# have an empty CHIP mutation match in the INFO field
# i.e. "CHIP_Transcript=.;" or "CHIP_Transcript=.,.;"
if [ "$(zgrep -v '^#' ${NEG_VCF} | awk '$7 ~ /chip_mutation_match_filter_fail|mutation_in_c_terminal/' | awk '$8 !~ /CHIP_Transcript=\.(,\.)*;/' | wc -l)" -ne 0 ]; then
    echo "FAIL: Negative sites VCF contains non-empty CHIP mutation matches in INFO field"
    exit 1
fi

# Test that all sites that have an empty CHIP mutation match in the INFO field
# are filtered out due to a CHIP mutation mismatch or not being a valid CHIP candidate
if [ "$(zgrep -v '^#' ${NEG_VCF} | awk '$8 ~ /CHIP_Transcript=\.(,\.)*;/' | awk '$7 !~ /chip_mutation_match_filter_fail|mutation_in_c_terminal|not_exonic_or_splicing_variant|exonic_or_splicing_variant_not_in_chip_transcript/' | wc -l)" -ne 0 ]; then
    echo "FAIL: Negative sites VCF contains non-CHIP mutation mismatch filter failures"
    exit 1
fi

# Test that all sites that have a non-empty CHIP mutation match in the INFO field
# are not filtered out due to a CHIP mutation mismatch or not being a valid CHIP candidate
if [ "$(zgrep -v '^#' ${NEG_VCF} | awk '$8 !~ /CHIP_Transcript=\.(,\.)*;/' | awk '$7 ~ /chip_mutation_match_filter_fail|mutation_in_c_terminal|not_exonic_or_splicing_variant|exonic_or_splicing_variant_not_in_chip_transcript/' | wc -l)" -ne 0 ]; then
    echo "FAIL: Negative sites VCF contains CHIP mutation mismatch filter failures"
    exit 1
fi

# === All-sites tests ===
# Check that the PASS sites match the positive sites
if [ "$(diff <(zgrep -v '^#' ${ALL_VCF} | awk '$7 == "PASS"' | sort) <(zgrep -v '^#' ${POS_VCF} | sort) | wc -l)" -ne 0 ]; then
    echo "FAIL: All sites VCF PASS sites do not match positive sites VCF"
    exit 1
fi

# Check that the failed sites match the negative sites
if [ "$(diff <(zgrep -v '^#' ${ALL_VCF} | awk '$7 != "PASS"' | sort) <(zgrep -v '^#' ${NEG_VCF} | sort) | wc -l)" -ne 0 ]; then
    echo "FAIL: All sites VCF failed sites do not match negative sites VCF"
    exit 1
fi

# === ALL TESTS PASSED ===
echo "All tests passed"
