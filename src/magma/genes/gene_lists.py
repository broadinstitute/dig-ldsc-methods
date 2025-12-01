import gzip
import json
from multiprocessing import Pool
import numpy as np
import os
import re
import subprocess
import time
from typing import Dict, List, Tuple
from zipfile import ZipFile

input_path = os.environ.get('INPUT_PATH')
s3_path = os.environ.get('S3_BUCKET')


def get_genes_path(data_path: str) -> str:
    return f'{data_path}/inputs/NCBI37.3.gene.loc'


def get_magma_path(data_path: str) -> str:
    return f'{data_path}/inputs/magma/magma'


def get_phenotype_genes_path(data_path: str, ancestry: str) -> str:
    return f'{data_path}/inputs/magma.{ancestry}.genes.zip'


def check_genes() -> None:
    if not os.path.exists(get_genes_path(input_path)):
        subprocess.check_call(f'./bootstrap/genes.bootstrap.sh {s3_path} {input_path}', shell=True)


def check_magma() -> None:
    if not os.path.exists(get_magma_path(input_path)):
        subprocess.check_call(f'./bootstrap/magma.bootstrap.sh {s3_path} {input_path}', shell=True)


def check_phenotype_genes(ancestry: str) -> None:
    if not os.path.exists(get_phenotype_genes_path(input_path, ancestry)):
        subprocess.check_call(f'./bootstrap/phenotype_genes.bootstrap.sh {s3_path} {input_path} {ancestry}', shell=True)


def get_gene_map() -> Dict:
    out = {}
    with open(f'{input_path}/inputs/NCBI37.3.gene.loc', 'r') as f:
        for line in f:
            split_line = [col.strip() for col in line.split('\t')]
            geneId = int(split_line[0])
            gene = split_line[5]
            out[gene] = geneId
    return out


def convert_gene_list(data_path: str, file: str) -> None:
    gene_map = get_gene_map()
    with open(f'{data_path}/raw/gene_list.txt', 'w') as f_out:
        with open(f'{data_path}/raw/{file}', 'r') as f:
            gene_ints = []
            for line in f:
                gene = line.strip()
                if gene in gene_map:
                    gene_ints.append(gene_map[gene])
            f_out.write('gene_list {}\n'.format(' '.join(map(str, gene_ints))))


def get_all_phenotypes(data_path: str, ancestry: str) -> List[str]:
    with ZipFile(get_phenotype_genes_path(data_path, ancestry)) as input_zip:
        return [phenotype for file in input_zip.filelist
                for phenotype in re.findall(r'([^.]*).genes.raw', file.filename)]


def run_phenotype(args: Tuple) -> List:
    data_path, phenotype, ancestry = args
    with ZipFile(get_phenotype_genes_path(input_path, ancestry)) as input_zip:
        t = time.time()
        input_zip.extract(f'{phenotype}.genes.raw', f'.')
        try:
            subprocess.check_call(f'{input_path}/inputs/magma/magma '
                                  f'--gene-results {phenotype}.genes.raw '
                                  f'--set-annot {data_path}/raw/gene_list.txt '
                                  f'--out {data_path}/magma/intermediate/{phenotype}.pathways', shell=True)
        except:
            return []
        os.remove(f'{phenotype}.genes.raw')
        print(phenotype, time.time() - t)
        return convert(data_path, phenotype)


def convert(data_path: str, phenotype: str) -> List:
    data = []
    with open(f'{data_path}/magma/intermediate/{phenotype}.pathways.gsa.out', 'r') as f:
        # Remove header
        line = f.readline().strip()
        while line[0] == '#':
            line = f.readline().strip()

        for line in f:
            split_line = [col.strip() for col in re.sub(' +', ' ', line).split(' ')]
            data.append({
                'phenotype': phenotype,
                'numGenes': int(split_line[2]),
                'beta': float(split_line[3]),
                'betaStdErr': float(split_line[4]),
                'stdErr': float(split_line[5]),
                'pValue': float(split_line[6]) if float(split_line[6]) > 0 else np.nextafter(0, 1)
            })
    return data


def save(data_path: str, data: List) -> None:
    with gzip.open(f'{data_path}/magma/gene_list/associations.gene_list.json.gz', 'wt') as f_out:
        for phenotype_data in data:
            for line_json in phenotype_data:
                f_out.write(json.dumps(line_json) + '\n')


def main(data_path: str, metadata: Dict) -> None:
    ancestry = metadata['ancestry']
    file = metadata['file']

    check_genes()
    check_magma()
    check_phenotype_genes(ancestry)

    convert_gene_list(data_path, file)

    os.makedirs(f'{data_path}/magma/intermediate', exist_ok=True)
    os.makedirs(f'{data_path}/magma/gene_list', exist_ok=True)
    inputs = [(data_path, phenotype, ancestry) for phenotype in get_all_phenotypes(input_path, ancestry)]
    with Pool(10) as p:
        data = p.map(run_phenotype, inputs)

    save(data_path, data)
