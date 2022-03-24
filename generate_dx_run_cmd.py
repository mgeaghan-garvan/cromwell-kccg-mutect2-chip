import json
import argparse
import re
import sys

def check_args(args=None):
    parser = argparse.ArgumentParser(description="Generate the run command for DNAnexus.", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-j", "--json", help="Path to input JSON file.", required=True, type=str)
    parser.add_argument("-b", "--batch", help="Path to a batch input TSV file for running multiple single-sample workflows in parallel. See the README and the template file at ./input/inputFiles.tsv for information on the required format.", default="", type=str)
    parser.add_argument("-d", "--destination", help="Path to output directory on the DNAnexus platform.", required=True, type=str)
    parser.add_argument("-w", "--workflow", help="Path to workflow on the DNAnexus platform.", required=True, type=str)
    parser.add_argument("-o", "--output", help="Path to bash file to output run command to. Default: _dx_run.sh.", default="_dx_run.sh", type=str)
    arguments = parser.parse_args(args)
    return arguments


VALID_JSON_KEYS = [
    "Mutect2CHIP",
    "Mutect2CHIP_CHIP",
    "Mutect2CHIP_Annovar",
    "Mutect2CHIP_VEP"
]

JSON_KEY_MAP = {
    "Mutect2CHIP": [
        "tumor_reads",
        "tumor_reads_index",
        "normal_reads",
        "normal_reads_index"
    ],
    "Mutect2CHIP_CHIP": [
        "input_vcf",
        "tumor_sample_name"
    ],
    "Mutect2CHIP_Annovar": [
        "input_vcf"
    ],
    "Mutect2CHIP_VEP": [
        "input_vcf",
        "input_vcf_idx"
    ]
}


def validate_batch_json(json_data):
    pass


def main(args):
    # Load JSON data
    with open(args.json, 'r') as f:
        data = json.load(f)

    # Replace the workflow name with "stage-common"
    new_data = {}
    for k in data.keys():
        new_key = re.sub("Mutect2CHIP(_[^\.]+)?\.", "stage-common.", k)
        new_data[new_key] = data[k]

    with open(args.output, 'w') as f:
        f.write("#!/bin/bash\n")
        f.write(f"dx run --watch --yes --destination {args.destination} \\\n    ")
        for k in new_data.keys():
            f.write(f"-i{k}=\"{new_data[k]}\" \\\n    ")
        f.write(args.workflow)
        f.write("\n")


if __name__ == "__main__":
    args = check_args(sys.argv[1:])
    main(args)
