import argparse
import gzip
import json
import os
import shutil
import subprocess
from typing import Dict, List

s3_bucket = 's3://dig-ldsc-server'
data_path = '.'


def download(username: str, dataset: str) -> Dict:
    path = f'{s3_bucket}/userdata/{username}/annotation/{dataset}/raw/'
    subprocess.check_call(f'aws s3 cp {path} dataset/ --recursive', shell=True)
    with open('dataset/metadata', 'r') as f:
        metadata = json.load(f)
    return metadata


def get_g1000_data(ancestry: str, chromosome: int) -> List:
    data = []
    with open(f'{data_path}/g1000/{ancestry}/chr{chromosome}.bim', 'r') as f:
        for line in f:
            data.append(int(line.strip().split('\t')[3]))
    return data


def get_file_data(file: str, chromosome: int) -> List:
    data = []
    with open(f'dataset/{file}', 'r') as f:
        for line in f:
            range_chr, start, end = line.strip().split('\t')[:3]
            if range_chr == str(chromosome):
                data.append((int(start), int(end)))
    return data


def write_annot(path: str, dataset: str, data: List, g1000_data: List):
    with gzip.open(path, 'wt') as f_out:
        f_out.write(f'{dataset}\n')
        idx = 0
        curr_start, curr_end = data[idx]
        for curr_position in iter(g1000_data):
            while idx < len(data) and curr_end < curr_position:
                idx += 1
                if idx < len(data):
                    curr_start, curr_end = data[idx]
            if curr_start < curr_position <= curr_end:
                f_out.write('1\n')
            else:
                f_out.write('0\n')


def run_chromosome(filename: str, dataset: str, ancestry: str, chromosome: int):
    data = get_file_data(filename, chromosome)
    g1000_data = get_g1000_data(ancestry, chromosome)

    if not os.path.exists('annotation'):
        os.mkdir('annotation')
    write_annot(f'annotation/{dataset}.{ancestry}.{chromosome}.annot.gz', dataset, data, g1000_data)


def zip_up_data(dataset: str, ancestry: str):
    if not os.path.exists('zip'):
        os.mkdir('zip')
    subprocess.check_call(f'zip -j zip/annotation.{dataset}.{ancestry}.zip annotation/{dataset}.{ancestry}.*.annot.gz', shell=True)


def upload(username: str, dataset: str):
    path = f'{s3_bucket}/userdata/{username}/annotation/{dataset}/ldsc/annotation/'
    subprocess.check_call(f'aws s3 cp zip/ {path} --recursive', shell=True)


def clean_up():
    for directory in ['dataset', 'annotation', 'zip']:
        if os.path.exists(directory):
            shutil.rmtree(directory)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--username', default=None, required=True, type=str)
    parser.add_argument('--dataset', default=None, required=True, type=str)
    args = parser.parse_args()

    metadata = download(args.username, args.dataset)
    for ancestry in ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']:
        for chromosome in range(1, 23):
            run_chromosome(metadata['file'], args.dataset, ancestry, chromosome)
        zip_up_data(args.dataset, ancestry)
    upload(args.username, args.dataset)
    clean_up()


if __name__ == '__main__':
    main()
