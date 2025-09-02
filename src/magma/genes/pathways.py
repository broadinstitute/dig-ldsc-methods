import argparse
import gzip
import json
import numpy as np
import os
import re
import shutil
import subprocess
from typing import Dict

input_path = os.environ.get('INPUT_PATH')
s3_path = os.environ.get('S3_BUCKET')


def get_pathways_path(data_path: str) -> str:
    return f'{data_path}/inputs/pathwayGenes.txt'


def get_magma_path(data_path: str) -> str:
    return f'{data_path}/inputs/magma/magma'


def check_pathways() -> None:
    if not os.path.exists(get_pathways_path(input_path)):
        subprocess.check_call(f'./bootstrap/pathways.bootstrap.sh {s3_path} {input_path}', shell=True)


def check_magma() -> None:
    if not os.path.exists(get_magma_path(input_path)):
        subprocess.check_call(f'./bootstrap/magma.bootstrap.sh {s3_path} {input_path}', shell=True)

def convert(data_path: str) -> None:
    with gzip.open(f'{data_path}/magma/genes/associations.pathways.json.gz', 'wt') as f_out:
        with open(f'{data_path}/magma/genes/associations.pathways.gsa.out', 'r') as f:
            # Remove header
            line = f.readline().strip()
            while line[0] == '#':
                line = f.readline().strip()

            for line in f:
                split_line = [col.strip() for col in re.sub(' +', ' ', line).split(' ')]
                pathway_name = split_line[7]
                dict_out = {
                    'pathwayName': pathway_name,
                    'numGenes': int(split_line[2]),
                    'beta': float(split_line[3]),
                    'betaStdErr': float(split_line[4]),
                    'stdErr': float(split_line[5]),
                    'pValue': float(split_line[6]) if float(split_line[6]) > 0 else np.nextafter(0, 1)
                }
                f_out.write(json.dumps(dict_out) + '\n')


def main(data_path: str) -> None:
    check_pathways()
    check_magma()

    os.makedirs(f'{data_path}/magma/genes', exist_ok=True)
    subprocess.check_call(f'{input_path}/inputs/magma/magma '
                          f'--gene-results {data_path}/magma/genes/associations.genes.raw '
                          f'--set-annot {input_path}/inputs/pathwayGenes.txt '
                          f'--out {data_path}/magma/genes/associations.pathways', shell=True)
    convert(data_path)
