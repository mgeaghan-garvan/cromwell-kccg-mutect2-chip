"""
Move the CHIP annotation outputs to the analysis bucket
"""
import hailtop.batch as hb
from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, image_path

import re

VALID_PROJECTS = [
    'kccg-genomics-med',
]


_config = get_config()
DATASET = _config['workflow']['dataset']
JOB_NAME = _config['workflow']['name']
WORKFLOW_NAME = _config['workflow']['workflow_name']
RUN_ID = _config['workflow']['run_id']
OUTPUT_PREFIX = f'{JOB_NAME}/{WORKFLOW_NAME}/{RUN_ID}'
CROMWELL_OUTPUT_PREFIX = f'mutect2-chip/{OUTPUT_PREFIX}'
CROMWELL_OUTPUT_BUCKET = _config['storage'][DATASET]['default']
CROMWELL_OUTPUT_PATH = f'{CROMWELL_OUTPUT_BUCKET}/{CROMWELL_OUTPUT_PREFIX}'
OUTPUT_BUCKET = _config['storage'][DATASET]['analysis']
OUTPUT_PATH = f'{OUTPUT_BUCKET}/{CROMWELL_OUTPUT_PREFIX}'

# Check for valid buckets
if CROMWELL_OUTPUT_BUCKET not in [
    f'gs://cpg-{project}-main'
    for project in VALID_PROJECTS
]:
    raise ValueError(f'Invalid CROMWELL_OUTPUT_BUCKET: {CROMWELL_OUTPUT_BUCKET}')
if OUTPUT_BUCKET not in [
    f'gs://cpg-{project}-main-analysis'
    for project in VALID_PROJECTS
]:
    raise ValueError(f'Invalid OUTPUT_BUCKET: {OUTPUT_BUCKET}')

b = get_batch()

process_j = b.new_job('move-chip-output-files')
process_j.image(image_path('cpg_workflows'))

FILES_TO_MOVE_LIST = []

if _config['analyses'].get('chip', False):
    FILES_TO_MOVE_LIST.extend([
        f'{CROMWELL_OUTPUT_PATH}/**/*.all_variants*.csv',
        f'{CROMWELL_OUTPUT_PATH}/**/*.exonic_splicing_variants.csv',
        f'{CROMWELL_OUTPUT_PATH}/**/*.chip*.csv',
    ])
if _config['analyses'].get('annovar', False):
    FILES_TO_MOVE_LIST.extend([
        f'{CROMWELL_OUTPUT_PATH}/**/*.hg38_multianno.vcf',
        f'{CROMWELL_OUTPUT_PATH}/**/*.hg38_multianno.txt',
    ])

# Check for malformed paths
for p in FILES_TO_MOVE_LIST:
    if not re.match(r'^gs://[A-Za-z0-9\-_\/\.]+$', p):
        raise ValueError(f'Invalid CROMWELL_OUTPUT_PATH: {p}')
if not re.match(r'^gs://[A-Za-z0-9\-_\/\.]+$', OUTPUT_PATH):
    raise ValueError(f'Invalid OUTPUT_PATH: {OUTPUT_PATH}')

FILES_TO_MOVE = ' '.join(FILES_TO_MOVE_LIST)

process_j.command(f"gcloud -q auth activate-service-account --key-file=/gsa-key/key.json; gsutil mv {FILES_TO_MOVE} {OUTPUT_PATH}")

b.run(wait=False)