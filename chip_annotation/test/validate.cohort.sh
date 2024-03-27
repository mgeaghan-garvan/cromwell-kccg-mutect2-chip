#!/bin/bash
# Validate test outputs

OUTDIR="output"

set -euo pipefail

# Check that the cohort filter ouptuts match the expected outputs
if [ "$(diff <(zgrep -h -v "^#" ${OUTDIR}/sample*.vcf.gz | sort -u ) test/data/sample1-5.outputs.tsv | wc -l)" -ne 0 ]; then
    echo "FAIL: All sites VCF failed sites do not match negative sites VCF"
    exit 1
fi

# === ALL TESTS PASSED ===
echo "All tests passed"
