import argparse
import gzip
import json
import numpy as np
import os
from scipy.stats import chi2
import shutil
import subprocess
from typing import List, Dict

var_id_columns = ['chromosome', 'position', 'reference', 'alt']
zn_columns = ['pValue', 'beta', 'n']

s3_bucket = 's3://dig-ldsc-server'
data_path = '.'


def download(username: str, dataset: str) -> Dict:
    path = f'{s3_bucket}/upload/{username}/{dataset}/genetic/'
    subprocess.check_call(f'aws s3 cp {path} dataset/ --recursive', shell=True)
    with open('dataset/metadata', 'r') as f:
        metadata = json.load(f)
    return metadata


def get_var_to_rs_map(ancestry: str) -> Dict:
    var_to_rs = {}
    with open(f'{data_path}/snplist/ldsc.{ancestry}.snplist', 'r') as f:
        for full_row in f.readlines():
            var_id, rs_id = full_row.strip().split('\t')
            var_to_rs[var_id] = rs_id
    return var_to_rs


def p_to_z(p: float, beta: float) -> float:
    return np.sqrt(chi2.isf(p, 1)) * (-1)**(beta < 0)


def valid_line(line: Dict, col_map: Dict) -> bool:
    return all([line.get(col_map[column]) is not None for column in var_id_columns]) and \
           all([line.get(col_map[column]) is not None for column in var_id_columns]) and \
           0 < float(line[col_map['pValue']]) <= 1


def stream_to_data(file: str, ancestry: str, col_map: Dict) -> List:
    var_to_rs_map = get_var_to_rs_map(ancestry)
    out = []
    count = 0
    with gzip.open(f'dataset/{file}', 'rt') as f_in:
        header = f_in.readline().strip().split('\t')
        for json_string in f_in:
            line = dict(zip(header, json_string.strip().split('\t')))
            count += 1
            if valid_line(line, col_map):
                chromosome, position, reference, alt = (line[col_map[c]] for c in var_id_columns)
                var_id = f'{chromosome}:{position}:{reference.upper()}:{alt.upper()}'
                if var_id in var_to_rs_map:
                    pValue, beta, n = (line[col_map[c]] for c in zn_columns)
                    out.append((var_to_rs_map[var_id], p_to_z(float(pValue), float(beta)), float(n)))
    print('{} SNPs translated from a total of {}'.format(len(out), count))
    return out


def filter_data_to_dict(data: List) -> Dict:
    N90 = np.quantile([a[2] for a in data], 0.9)
    return {d[0]: (d[1], d[2]) for d in data if d[2] >= N90 / 1.5}


def save_to_file(dataset: str, ancestry: str, data: Dict) -> str:
    if not os.path.exists('sumstats'):
        os.mkdir('sumstats')
    out_file = f'sumstats/{dataset}.sumstats.gz'
    with gzip.open(out_file, 'wt') as f:
        f.write('SNP\tZ\tN\n')
        for rs_id in ld_rs_iter(ancestry):
            if rs_id in data:
                f.write('{}\t{}\t{}\n'.format(rs_id, round(data[rs_id][0], 3), data[rs_id][1]))
            else:
                f.write(f'{rs_id}\t\t\n')
    return out_file


def ld_rs_iter(ancestry: str) -> str:
    for CHR in range(1, 23):
        with gzip.open(f'{data_path}/weights/{ancestry}/weights.{CHR}.l2.ldscore.gz', 'rt') as f:
            _ = f.readline()
            for line in f:
                yield line.strip().split('\t')[1]


def upload(username: str, dataset: str, file: str):
    path = f'{s3_bucket}/data/{username}/{dataset}/ldsc/sumstats/'
    subprocess.check_call(f'aws s3 cp {file} {path}', shell=True)


def clean_up():
    for directory in ['dataset', 'sumstats']:
        if os.path.exists(directory):
            shutil.rmtree(directory)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--username', default=None, required=True, type=str)
    parser.add_argument('--dataset', default=None, required=True, type=str)
    args = parser.parse_args()

    metadata = download(args.username, args.dataset)
    data = stream_to_data(metadata['file'], metadata['ancestry'], metadata['col_map'])
    if len(data) > 0:
        data_dict = filter_data_to_dict(data)
        out_file = save_to_file(args.dataset, metadata['ancestry'], data_dict)
        upload(args.username, args.dataset, out_file)
    clean_up()


if __name__ == '__main__':
    main()
