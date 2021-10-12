import json
import subprocess

with open("run_id.txt", 'r') as f:
    id = json.load(f).get('id', None)

if id is None:
    raise ValueError("Workflow ID is not present in run_id.txt")

subprocess.run(['./abort.sh', id])
