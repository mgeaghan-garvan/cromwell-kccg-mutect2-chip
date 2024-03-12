import json
import argparse
import re
import sys
import csv
import os
import pandas as pd

def check_args(args=None):
    parser = argparse.ArgumentParser(description="Generate batch input JSON files from a single completed JSON input file and a TSV file detailing the variable parameters for each sample.", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-j", "--json", help="Path to input JSON file.", required=True, type=str)
    parser.add_argument("-b", "--batch", help="Path to a batch input TSV file. Examples can be found within input/inputFiles/batch", default="", type=str)
    parser.add_argument("-o", "--output", help="Directory into which the batch input JSON files should be placed.", default="input/config/batch", type=str)
    parser.add_argument("-f", "--format", help="Output format, one of either 'singlefile' (default), 'multifile', or 'dx'. 'singlefile' will output a single JSON file containing an array of JSON objects, one per sample. 'multifile' will output multiple JSON files, one per sample. 'dx' will output multiple DNAnexus-ready input files, one per sample.", default='singlefile', choices=['singlefile', 'multifile', 'dx'], type=str)
    parser.add_argument("-d", "--dx_destination", help="DNAnexus output path. Required when using -f/--format dx.", default="", type=str)
    parser.add_argument("-p", "--dx_priority", help="Priority for DNAnexus run. One of 'low' (default), 'normal', or 'high'.", default="low", type=str)
    parser.add_argument("-w", "--dx_workflow", help="DNAnexus workflow path. Required when using -f/--format dx.", default="", type=str)
    arguments = parser.parse_args(args)
    if arguments.format == 'dx' and (
        not arguments.dx_destination or
        not arguments.dx_workflow
    ):
        raise ValueError("ERROR: Must supply -d/--dx_destination and -w/--dx_workflow when using -f/--format dx.")
    return arguments


VALID_KEYS = {
    "Mutect2CHIP": [
        "tumor_reads",
        "tumor_reads_index",
        "normal_reads",
        "normal_reads_index",
    ],
    "Mutect2": [
        "tumor_reads",
        "tumor_reads_index",
        "normal_reads",
        "normal_reads_index",
    ],
    "CallU2AF1": [
        "tumor_reads",
        "tumor_reads_index",
        "mutect2_output_vcf",
        "mutect2_output_vcf_index",
    ],
    "CHIP": [
        "tumor_sample_name",
        "input_vcf",
    ],
    "Annovar": [
        "input_vcf",
    ],
    "SpliceAI": [
        "input_vcf",
    ],
    "VEP": [
        "input_vcf",
        "input_vcf_idx",
    ],
}


OPTIONAL_KEYS = {
    "Mutect2CHIP": [
        "normal_reads",
        "normal_reads_index",
    ],
    "Mutect2": [
        "normal_reads",
        "normal_reads_index",
    ],
}


def validate_keys(json_data, batch_data):
    json_keys = list(json_data.keys())
    json_keys_split = [k.split(".") for k in json_keys]
    json_keys_prefix = [k[0] for k in json_keys_split]
    json_keys_prefix = list(set(json_keys_prefix))
    if len(json_keys_prefix) != 1 or json_keys_prefix[0] not in VALID_KEYS:
        raise ValueError("ERROR: Invalid input JSON file.")
    json_key_prefix = json_keys_prefix[0]
    all_valid_keys = VALID_KEYS.get(json_key_prefix)
    all_optional_keys = OPTIONAL_KEYS.get(json_key_prefix, [])
    all_required_keys = [k for k in all_valid_keys if k not in all_optional_keys]
    batch_header = list(batch_data.columns)
    if any([h not in all_valid_keys for h in batch_header if h != 'sample_id']):
        raise ValueError(f"ERROR: Invalid column headers defined in TSV file: {batch_header}")
    if any([k not in batch_header for k in all_required_keys]):
        raise ValueError(f"ERROR: Required column headers missing from TSV file: {all_required_keys}")
    batch_parameters = []
    data = batch_data.fillna('').transpose().to_dict()
    for sample_id in data:
        sample = data[sample_id]
        sample_parameters = {'sample_id': sample['sample_id'], 'params': {}}
        for k in all_valid_keys:
            json_key = f"{json_key_prefix}.{k}"
            optional = k in all_optional_keys
            value = sample.get(k, '')
            if not value and not optional:
                raise ValueError(f"ERROR: Required value missing from TSV file: {all_required_keys}")
            sample_parameters['params'][json_key] = {
                'value': value,
                'optional': optional,
            }
        batch_parameters.append(sample_parameters)
    return batch_parameters


def convert_to_dx(json_data):
    new_json_data = {f"stage-common.{k.split('.')[1]}": v for k, v in json_data.items()}
    return new_json_data


def main(args):
    # Load JSON data
    with open(args.json, 'r') as f:
        config = json.load(f)
    json_basename = os.path.splitext(os.path.basename(args.json))[0]
    if not args.batch and args.format == 'dx':
        # DNAnexus jobs aren't required to be batched, so if no batch file is provided,
        # just convert the single input JSON to a DNAnexus input file and exit
        sample_json = config.copy()
        sample_dx = convert_to_dx(sample_json)
        dx_out = os.path.join(args.output, f"{json_basename}.dx.sh")
        if os.path.exists(dx_out):
            raise FileExistsError(f"ERROR: DNAnexus submission script file {dx_out} already exists!")
        with open(dx_out, 'w') as f:
            f.write("#!/bin/bash\n")
            dx_destination = re.sub("/$", "", args.dx_destination)
            f.write(f"dx run --priority {args.dx_priority} --yes --destination {dx_destination}/ \\\n    ")
            for k, v in sample_dx.items():
                f.write(f"-i{k}=\"{v}\" \\\n    ")
            f.write(args.dx_workflow)
            f.write("\n")
        return
    elif not args.batch:
        # All other formats require a batch file, so if none is provided, exit
        raise ValueError("ERROR: Must supply -b/--batch when using -f/--format singlefile or multifile.")
    batch_data = pd.read_table(args.batch, sep = '\s+')
    batch_params = validate_keys(config, batch_data)
    singlefile_array=[]
    for sample in batch_params:
        sample_id = sample['sample_id']
        sample_json = config.copy()
        for param in sample['params']:
            param_value = sample['params'][param]['value']
            param_optional = sample['params'][param]['optional']
            if not param_value and param_optional:
                sample_json.pop(param, None)
                continue
            sample_json[param] = param_value
        if args.format == 'multifile':
            out_file = os.path.join(args.output, f"{json_basename}.{sample_id}.json")
            if os.path.exists(out_file):
                raise FileExistsError(f"ERROR: JSON file {out_file} already exists!")
            with open(out_file, 'w') as f:
                json.dump(sample_json, f)
        elif args.format == 'singlefile':
            singlefile_array.append(sample_json)
        elif args.format == 'dx':
            sample_dx = convert_to_dx(sample_json)
            dx_out = os.path.join(args.output, f"{json_basename}.{sample_id}.dx.sh")
            if os.path.exists(dx_out):
                raise FileExistsError(f"ERROR: DNAnexus submission script file {dx_out} already exists!")
            with open(dx_out, 'w') as f:
                f.write("#!/bin/bash\n")
                dx_destination = re.sub("/$", "", args.dx_destination)
                f.write(f"dx run --priority {args.dx_priority} --yes --destination {dx_destination}/{sample_id}/ \\\n    ")
                for k, v in sample_dx.items():
                    f.write(f"-i{k}=\"{v}\" \\\n    ")
                f.write(args.dx_workflow)
                f.write("\n")
    if singlefile_array:
        singlefile_out = os.path.join(args.output, f"{json_basename}.batch.json")
        if os.path.exists(singlefile_out):
            raise FileExistsError(f"ERROR: JSON file {singlefile_out} already exists!")
        with open(singlefile_out, 'w') as f:
            json.dump(singlefile_array, f)


if __name__ == "__main__":
    args = check_args(sys.argv[1:])
    main(args)
