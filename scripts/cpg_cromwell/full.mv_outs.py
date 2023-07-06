"""
Move the analysis outputs of the full pipeline to the analysis bucket
"""
import hailtop.batch as hb
from cpg_utils.config import get_config
from cpg_utils.hail_batch import output_path, remote_tmpdir
from analysis_runner.cromwell import (
    run_cromwell_workflow_from_repo_and_get_outputs,
    CromwellOutputType,
)
import toml


_config = get_config()
BILLING_PROJECT = _config['hail']['billing_project']
DATASET = _config['workflow']['dataset']
ACCESS_LEVEL = _config['workflow']['access_level']
JOB_NAME = _config['workflow']['name']
WORKFLOW_NAME = _config['workflow']['workflow_name']
RUN_ID = _config['workflow']['run_id']
OUTPUT_PREFIX = f'{JOB_NAME}/{WORKFLOW_NAME}/{RUN_ID}'
CROMWELL_OUTPUT_PREFIX = f'mgeaghan/mutect2-chip/{OUTPUT_PREFIX}'
CROMWELL_OUTPUT_BUCKET = _config['storage'][DATASET]['default']
CROMWELL_OUTPUT_PATH = f'{CROMWELL_OUTPUT_BUCKET}/{CROMWELL_OUTPUT_PREFIX}'
OUTPUT_BUCKET = _config['storage'][DATASET]['analysis']
OUTPUT_PATH = f'{OUTPUT_BUCKET}/{CROMWELL_OUTPUT_PREFIX}'

sb = hb.ServiceBackend(billing_project=BILLING_PROJECT, remote_tmpdir=remote_tmpdir())
b = hb.Batch(backend=sb, default_image=_config['workflow']['driver_image'])

process_j = b.new_job('move-cromwell-output-files')

# TODO: in progress
FILES_TO_MOVE_LIST = []
if _config['move_outputs'].get('vcf', False):
    FILES_TO_MOVE_LIST.append(f'{CROMWELL_OUTPUT_PATH}/**/*.vcf')
    FILES_TO_MOVE_LIST.append(f'{CROMWELL_OUTPUT_PATH}/**/*.vcf.gz')
    FILES_TO_MOVE_LIST.append(f'{CROMWELL_OUTPUT_PATH}/**/*.vcf.gz.tbi')
if _config['move_outputs'].get('bam', False):
    FILES_TO_MOVE_LIST.append(f'{CROMWELL_OUTPUT_PATH}/**/*.bam')
    FILES_TO_MOVE_LIST.append(f'{CROMWELL_OUTPUT_PATH}/**/*.bai')
if _config['move_outputs'].get('annovar', False):
    FILES_TO_MOVE_LIST.append(f'{CROMWELL_OUTPUT_PATH}/**/*_multianno.txt')
if _config['move_outputs'].get('chip', False):
    FILES_TO_MOVE_LIST.append(f'{CROMWELL_OUTPUT_PATH}/**/*.all_variants*.csv')
    FILES_TO_MOVE_LIST.append(f'{CROMWELL_OUTPUT_PATH}/**/*.exonic_splicing_variants.csv')
    FILES_TO_MOVE_LIST.append(f'{CROMWELL_OUTPUT_PATH}/**/*.chip*.csv')
FILES_TO_MOVE = ' '.join(FILES_TO_MOVE_LIST)

process_j.command(f"gcloud -q auth activate-service-account --key-file=/gsa-key/key.json; gsutil mv {FILES_TO_MOVE} {OUTPUT_PATH}")

b.run(wait=False)