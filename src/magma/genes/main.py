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


def check_envvars():
    assert input_path is not None
    assert s3_path is not None


def get_genes_annot_path(data_path: str) -> str:
    return f'{data_path}/inputs/all.genes.annot'


def get_genes_path(data_path: str) -> str:
    return f'{data_path}/inputs/NCBI37.3.gene.loc'


def get_magma_path(data_path: str) -> str:
    return f'{data_path}/inputs/magma/magma'


def get_g1000_path(data_path: str, ancestry: str) -> str:
    return f'{data_path}/inputs/g1000/{ancestry}/'


def check_genes_annot() -> None:
    if not os.path.exists(get_genes_annot_path(input_path)):
        subprocess.check_call(f'./bootstrap/genes_annot.bootstrap.sh {s3_path} {input_path}', shell=True)


def check_genes() -> None:
    if not os.path.exists(get_genes_path(input_path)):
        subprocess.check_call(f'./bootstrap/genes.bootstrap.sh {s3_path} {input_path}', shell=True)


def check_magma() -> None:
    if not os.path.exists(get_magma_path(input_path)):
        subprocess.check_call(f'./bootstrap/magma.bootstrap.sh {s3_path} {input_path}', shell=True)


def check_g1000(ancestry: str) -> None:
    if not os.path.exists(get_g1000_path(input_path, ancestry)):
        subprocess.check_call(f'./bootstrap/g1000.bootstrap.sh {s3_path} {input_path} {ancestry}', shell=True)


def get_metadata(data_path: str) -> Dict:
    with open(f'{data_path}/sumstats/metadata', 'r') as f:
        metadata = json.load(f)
    return metadata


def unzip_sumstats(data_path: str) -> None:
    with gzip.open(f'{data_path}/sumstats/magma.sumstats.csv.gz', 'rb') as f_in:
        with open(f'{data_path}/sumstats/magma.sumstats.csv', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)


def get_gene_map() -> Dict:
    out = {}
    with open(f'{input_path}/inputs/NCBI37.3.gene.loc', 'r') as f:
        for line in f:
            split_line = [col.strip() for col in line.split('\t')]
            geneId = int(split_line[0])
            gene = split_line[5]
            out[geneId] = gene
    return out


def convert(data_path: str) -> None:
    gene_map = get_gene_map()
    with gzip.open(f'{data_path}/genes/associations.genes.json.gz', 'wt') as f_out:
        with open(f'{data_path}/genes/associations.genes.out', 'r') as f:
            _ = f.readline()
            for line in f:
                split_line = re.sub(' +', ' ', line).split(' ')
                dict_out = {
                    'gene': gene_map[int(split_line[0])],
                    'nParam': int(split_line[5]),
                    'subjects': int(split_line[6]),
                    'zStat': float(split_line[7]),
                    'pValue': float(split_line[8]) if float(split_line[8]) > 0 else np.nextafter(0, 1)
                }
                f_out.write(json.dumps(dict_out) + '\n')


def main() -> None:
    check_envvars()
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir', default=None, required=True, type=str)
    data_path = parser.parse_args().dir

    metadata = get_metadata(data_path)
    ancestry = metadata['ancestry']

    check_genes_annot()
    check_genes()
    check_magma()
    check_g1000(ancestry)

    unzip_sumstats(data_path)
    os.makedirs(f'{data_path}/genes', exist_ok=True)
    subprocess.check_call(f'{input_path}/inputs/magma/magma '
                          f'--bfile {input_path}/inputs/g1000/EUR/g1000_{ancestry} '
                          f'--pval {data_path}/sumstats/magma.sumstats.csv ncol=N '
                          f'--gene-annot {input_path}/inputs/all.genes.annot '
                          f'--out {data_path}/genes/associations', shell=True)
    convert(data_path)


if __name__ == '__main__':
    main()