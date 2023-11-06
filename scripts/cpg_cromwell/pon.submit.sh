analysis-runner \
  --dataset kccg-genomics-med \
  --description "PoN test run" \
  --output-dir "mgeaghan/submit/pon" \
  --access-level test \
  --config /Users/micgea/Documents/projects/chip/cohorts/tob/cpg_gcp_analysis/231106_wgs_full_1049_reanalysis/cromwell-kccg-mutect2-chip/input/inputs.pon.toml \
  --image australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_workflows:latest \
  python3 scripts/cpg_cromwell/pon.submit.py