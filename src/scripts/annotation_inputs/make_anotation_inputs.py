import glob
import gzip
import json
import numpy as np
from numpy import typing as npt
import os
import re
from typing import List

import weights, xtx_xty, sumstats


def baseline_annotation_path(ancestry: str, chromosome: int) -> str:
    return f'baseline/{ancestry}/baselineLD.{chromosome}.annot.gz'

def baseline_path(ancestry: str, chromosome: int) -> str:
    return f'baseline/{ancestry}/baselineLD.{chromosome}.l2.ldscore.gz'

def baseline_m_path(ancestry: str, chromosome: int) -> str:
    return f'baseline/{ancestry}/baselineLD.{chromosome}.l2.M_5_50'

def frq_path(ancestry: str, chromosome: int) -> str:
    return f'frq/{ancestry}/chr.{chromosome}.frq'

def get_baseline_ld_score(ancestry: str) -> npt.NDArray:
    output = []
    for chromosome in range(1, 23):
        with gzip.open(baseline_path(ancestry, chromosome), 'rt') as f:
            _ = f.readline()  # unused header
            for line in f:
                output.append(list(map(float, line.strip().split('\t')[3:])))
    return np.array(output)

def get_baseline_variables(ancestry: str) -> List[str]:
    with gzip.open(baseline_path(ancestry, 1), 'rt') as f:
        return [f'baseline___{column}' for column in f.readline().strip().split('\t')[3:]]

def get_baseline_parameter_snps(ancestry: str) -> npt.NDArray:
    output = []
    for chromosome in range(1, 23):
        with open(baseline_m_path(ancestry, chromosome), 'r') as f:
            output.append(list(map(float, f.readline().strip().split())))
    return np.array(output).sum(axis=0, keepdims=True).T

def get_baseline_annot(ancestry: str) -> npt.NDArray:
    out = []
    for chromosome in range(1, 23):
        with gzip.open(baseline_annotation_path(ancestry, chromosome), 'rt') as baseline_file, \
                open(frq_path(ancestry, chromosome), 'r') as frq_file:
            baseline_file.readline(), frq_file.readline()
            for baseline_line, frq_line in zip(baseline_file, frq_file):
                frq = float(frq_line.strip().split('\t')[4])
                if 0.05 < frq < 0.95:
                    values = baseline_line.strip().split('\t')[4:]
                    out.append(list(map(float, values)))
    return np.array(out)

def get_heritability(ancestry: str):
    heritability_z = {}
    for file in glob.glob('heritability/*/*.json'):
        with open(file, 'r') as f:
            for line in f:
                json_line = json.loads(line.strip())
                if json_line['ancestry'] == ancestry:
                    heritability_z[json_line['phenotype']] = float(json_line['h2'] / json_line['stdErr'])
    return heritability_z

def main():
    ancestry = 'EUR'
    heritability_z = get_heritability('EU')

    os.makedirs(f'inputs/{ancestry}/baseline', exist_ok=True)
    baseline_annot = get_baseline_annot(ancestry)
    baseline_ld = get_baseline_ld_score(ancestry)
    baseline_variables = get_baseline_variables(ancestry)
    baseline_parameter_snps = get_baseline_parameter_snps(ancestry)

    np.save(f'inputs/{ancestry}/baseline/baseline_annot.npy', baseline_annot)
    np.save(f'inputs/{ancestry}/baseline/baseline_ld.npy', baseline_ld)
    with open(f'inputs/{ancestry}/baseline/baseline_variables.txt', 'w') as f:
        f.write('\t'.join(baseline_variables))
    np.save(f'inputs/{ancestry}/baseline/baseline_parameter_snps.npy', baseline_parameter_snps)

    input_weights = weights.get_input_weights('.', ancestry)
    total_files = len(glob.glob(f'sumstats/{ancestry}/*.sumstats.gz'))
    for i, data_path in enumerate(glob.glob(f'sumstats/{ancestry}/*.sumstats.gz')):
        phenotype = re.findall(r'.*/([^/]*).sumstats.gz', data_path)[0].rsplit('_', 1)[0]
        print(phenotype, f'{i} / {total_files}', heritability_z.get(phenotype, 0.0))
        if heritability_z.get(phenotype, 0.0) > 7:
            os.makedirs(f'inputs/{ancestry}/phenotypes/{phenotype}', exist_ok=True)
            chisq, sample_size, idxs = sumstats.load_sumstats_from(data_path)
            chisq, sample_size, idxs = sumstats.filter_sumstats(chisq, sample_size, idxs)
            if idxs.shape[0] > 200000:
                sumstats_weights = input_weights[idxs]
                baseline_sumstats_ld = baseline_ld[idxs, :]

                weight = xtx_xty.get_weight(baseline_sumstats_ld, sumstats_weights, sample_size, chisq, baseline_parameter_snps)
                y = xtx_xty.get_y(chisq, weight)
                np.save(f'inputs/{ancestry}/phenotypes/{phenotype}/weights.npy', weight)
                np.save(f'inputs/{ancestry}/phenotypes/{phenotype}/sample_size.npy', sample_size)
                np.save(f'inputs/{ancestry}/phenotypes/{phenotype}/y.npy', y)
                np.save(f'inputs/{ancestry}/phenotypes/{phenotype}/idxs.npy', idxs)
            else:
                print(f'Skipping low SNP count phenotype {phenotype}')
        else:
            print(f'Skipping low heritability phenotype {phenotype}')

if __name__ == '__main__':
    main()
