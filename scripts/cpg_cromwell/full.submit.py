"""
Submit the full pipeline to cromwell via hail batch
"""
import hailtop.batch as hb
from hailtop.batch import Resource
from hailtop.batch.job import Job
from cpg_utils.config import get_config
from cpg_utils.hail_batch import remote_tmpdir
from analysis_runner.cromwell import (
    run_cromwell_workflow,
)
import toml
import os
from typing import List, Dict, Optional, Any
from analysis_runner.git import (
    get_git_default_remote,
    get_git_commit_ref_of_current_repository,
    get_repo_name_from_remote,
    prepare_git_job,
)


def submit_cromwell_workflow(
    b,
    job_prefix: str,
    dataset: str,
    access_level,
    workflow: str,
    libs: List[str],
    output_prefix: str,
    labels: Dict[str, str] = None,
    input_dict: Optional[Dict[str, Any]] = None,
    input_paths: List[str] = None,
    repo: Optional[str] = None,
    commit: Optional[str] = None,
    cwd: Optional[str] = None,
    driver_image: Optional[str] = None,
    project: Optional[str] = None,
    copy_outputs_to_gcp: bool = True,
) -> tuple[Job, dict[str, Resource]]:
    _driver_image = driver_image or os.getenv('DRIVER_IMAGE')

    submit_job = b.new_job(f'{job_prefix}_submit')
    submit_job.image(_driver_image)
    prepare_git_job(
        job=submit_job,
        repo_name=(repo or get_repo_name_from_remote(get_git_default_remote())),
        commit=(commit or get_git_commit_ref_of_current_repository()),
        is_test=access_level == 'test',
    )

    workflow_id_file = run_cromwell_workflow(
        job=submit_job,
        dataset=dataset,
        access_level=access_level,
        workflow=workflow,
        cwd=cwd,
        libs=libs,
        output_prefix=output_prefix,
        input_dict=input_dict,
        input_paths=input_paths,
        labels=labels,
        project=project,
        copy_outputs_to_gcp=copy_outputs_to_gcp,
    )

    return submit_job, workflow_id_file


_config = get_config()
BILLING_PROJECT = _config['hail']['billing_project']
DATASET = _config['workflow']['dataset']
ACCESS_LEVEL = _config['workflow']['access_level']
JOB_NAME = _config['workflow']['name']
OUTPUT_PREFIX = f'mutect2-chip/{JOB_NAME}'
DRIVER_IMAGE = _config['workflow']['driver_image']

sb = hb.ServiceBackend(billing_project=BILLING_PROJECT, remote_tmpdir=remote_tmpdir())
b = hb.Batch(backend=sb, default_image=DRIVER_IMAGE)

input_prefix = 'Mutect2CHIP'

input_dict = {
    f'{input_prefix}.{k}': v
    for k, v in _config['mutect2_chip']
}

submit_j, workflow_id_file = submit_cromwell_workflow(
    b=b,
    job_prefix='mutect2-chip-full',
    dataset=DATASET,
    access_level=ACCESS_LEVEL,
    workflow='full.wdl',
    output_prefix=OUTPUT_PREFIX,
    cwd='workflow/',
    input_dict=input_dict,
    libs=[],
    copy_outputs_to_gcp=True,
    driver_image=DRIVER_IMAGE,
)

b.run(wait=False)