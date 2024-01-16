"""
Remove the CHIP pipeline outputs from the release bucket once no longer needed
"""
import hailtop.batch as hb
from cpg_utils.config import get_config
from cpg_utils.hail_batch import remote_tmpdir, get_batch

import re

_config = get_config()
DATASET = _config['workflow']['dataset']
JOB_NAME = _config['workflow']['name']
WORKFLOW_NAME = _config['workflow']['workflow_name']
RUN_ID = _config['workflow']['run_id']
OUTPUT_PREFIX = f'{JOB_NAME}/{WORKFLOW_NAME}/{RUN_ID}'
CROMWELL_OUTPUT_PREFIX = f'mutect2-chip/{OUTPUT_PREFIX}'
CROMWELL_OUTPUT_BUCKET = _config['storage'][DATASET]['default']
# OUTPUT_BUCKET = _config['storage'][DATASET]['release']  # Currently no 'release' key in config TOML
OUTPUT_BUCKET = re.sub(r'\-main$', '-release', CROMWELL_OUTPUT_BUCKET)
OUTPUT_PATH = f'{OUTPUT_BUCKET}/{CROMWELL_OUTPUT_PREFIX}'

b = get_batch()

process_j = b.new_job('move-chip-output-files')

process_j.command(f"gcloud -q auth activate-service-account --key-file=/gsa-key/key.json; gsutil rm -r {OUTPUT_PATH}")

b.run(wait=False)