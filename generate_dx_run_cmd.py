import json
import argparse
import re
import sys
import csv
import os

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
    json_keys = list(json_data.keys())
    json_keys = [k.split(".")[0] for k in json_keys]
    json_keys = list(set(json_keys))
    if len(json_keys) != 1 or json_keys[0] not in VALID_JSON_KEYS:
        raise ValueError("ERROR: Invalid input JSON file.")
    return JSON_KEY_MAP[json_keys[0]]


def main(args):
    # Load JSON data
    with open(args.json, 'r') as f:
        data = json.load(f)

    json_data_list = []
    out_file_list = []
    if args.batch:
        batch_json_keys = validate_batch_json(data)
        with open(args.batch, 'r') as f:
            reader = csv.reader(f, delimiter="\t")
            i = 1
            for line in reader:
                if len(line) > len(batch_json_keys):
                    raise ValueError("ERROR: Invalid batch input TSV file.")
                new_json_data = data.copy()
                for idx, val in enumerate(line):
                    new_json_data[batch_json_keys[idx]] = val
                json_data_list.append(new_json_data)
                out_file_split = os.path.split(args.output)
                out_file = f"{out_file_split[0]}/_batch_{str(i)}_{out_file_split[1]}"
                out_file_list.append(out_file)
    else:
        json_data_list = [data]
        out_file_list = [args.output]

    for i, j in enumerate(json_data_list):
        # Replace the workflow name with "stage-common"
        new_data = {}
        for k in j.keys():
            new_key = re.sub("Mutect2CHIP(_[^\.]+)?\.", "stage-common.", k)
            new_data[new_key] = j[k]

        out_file = out_file_list[i]
        with open(out_file, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write(f"dx run --watch --yes --destination {args.destination} \\\n    ")
            for k in new_data.keys():
                f.write(f"-i{k}=\"{new_data[k]}\" \\\n    ")
            f.write(args.workflow)
            f.write("\n")


if __name__ == "__main__":
    args = check_args(sys.argv[1:])
    main(args)
