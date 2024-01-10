"""
Move the CHIP pipeline outputs to the tmp bucket for export to the GML private buckets
"""
import hailtop.batch as hb
from cpg_utils.config import get_config
from cpg_utils.hail_batch import remote_tmpdir, get_batch


_config = get_config()
DATASET = _config['workflow']['dataset']
JOB_NAME = _config['workflow']['name']
WORKFLOW_NAME = _config['workflow']['workflow_name']
RUN_ID = _config['workflow']['run_id']
OUTPUT_PREFIX = f'{JOB_NAME}/{WORKFLOW_NAME}/{RUN_ID}'
CROMWELL_OUTPUT_PREFIX = f'mutect2-chip/{OUTPUT_PREFIX}'
CROMWELL_OUTPUT_BUCKET = _config['storage'][DATASET]['default']
CROMWELL_OUTPUT_PATH = f'{CROMWELL_OUTPUT_BUCKET}/{CROMWELL_OUTPUT_PREFIX}'
OUTPUT_BUCKET = _config['storage'][DATASET]['tmp']
OUTPUT_PATH = f'{OUTPUT_BUCKET}/{CROMWELL_OUTPUT_PREFIX}'

b = get_batch()

process_j = b.new_job('move-chip-output-files')

process_j.command(f"gcloud -q auth activate-service-account --key-file=/gsa-key/key.json; gsutil mv {CROMWELL_OUTPUT_PATH} {OUTPUT_PATH}")

b.run(wait=False)