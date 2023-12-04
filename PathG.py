"""
Create by Tianyng W	2023/11/28
"""

import os
import json
import argparse
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(description='''The help of this script
======================================================================
Usage:
    python3 PathG.py -d <path> -o <output>
This script is used to generate the json file of sample fasq files list
these file should end with `.fq.gz` or `.fastq.gz`
Cretae by Tianyng W	2023/11/28
Bug report:https://github.com/heihaheihaha/Call_variants_V1.0.0/issues
======================================================================''', formatter_class=RawTextHelpFormatter)
parser.add_argument('--directory', '-d', required=True, help='The directory containing the fastq files')
parser.add_argument('--output', '-o', required=False, help='The Path of output JSON file, default: All_sample.json')

args = parser.parse_args()


import os
import json

def extract_details(filename):
    """
    Extracts sample name, lane, and read number from a given fastq file name.
    Example: 'V300102848_L03_86_1.fq.gz' -> ('V300102848', 'L03', '1')
    """
    parts = filename.split('_')
    sample_name = parts[0]
    lane = parts[1]
    read = filename.split('_')[-1].split('.')[0]
    return sample_name, lane, read

def process_directory(directory):
    """
    Traverses the directory and organizes fastq files into a JSON structure grouped by sample names.
    """
    samples = {}
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith('.fq.gz') or file.endswith('.fastq.gz'):
                sample_name, lane, read = extract_details(file)
                if sample_name not in samples:
                    samples[sample_name] = {}
                if lane not in samples[sample_name]:
                    samples[sample_name][lane] = {'R1': '', 'R2': ''}
                read_key = 'R1' if read == '1' else 'R2'
                samples[sample_name][lane][read_key] = os.path.join(root, file)
    return samples

if __name__ == '__main__':

	directory_path = args.directory
	sample_data = process_directory(directory_path)
	json_data = json.dumps(sample_data, indent=4)
	print(json_data)

	if not args.output:
		args.output = 'All_sample.json'

	sample_data = process_directory(args.directory)
	json_data = json.dumps(sample_data, indent=4)
	with open(args.output, 'w') as f:
		f.write(json_data)