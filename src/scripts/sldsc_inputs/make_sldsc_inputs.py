import glob
import gzip
import numpy as np
from numpy import typing as npt
import os
import re
from typing import List


def baseline_path(ancestry: str, chromosome: int) -> str:
    return f'baseline/{ancestry}/baselineLD.{chromosome}.l2.ldscore.gz'


def baseline_m_path(ancestry: str, chromosome: int) -> str:
    return f'baseline/{ancestry}/baselineLD.{chromosome}.l2.M_5_50'


def tissue_path(tissue: str, ancestry: str, chromosome: int) -> str:
    return f'tissue/{ancestry}/{tissue}/tissueLD.{chromosome}.l2.ldscore.gz'


def tissue_m_path(tissue: str, ancestry: str, chromosome: int) -> str:
    return f'tissue/{ancestry}/{tissue}/tissueLD.{chromosome}.l2.M_5_50'


def get_baseline_variables(ancestry: str) -> List[str]:
    with gzip.open(baseline_path(ancestry, 1), 'rt') as f:
        return [f'baseline___{column}' for column in f.readline().strip().split('\t')[3:]]


def get_tissue_variables(tissue: str, ancestry: str) -> List[str]:
    with gzip.open(tissue_path(tissue, ancestry, 1), 'rt') as f:
        return [f'tissue___{column}' for column in f.readline().strip().split('\t')[3:]]


def get_baseline_parameter_snps(ancestry: str) -> List[str]:
    output = []
    for chromosome in range(1, 23):
        with open(baseline_m_path(ancestry, chromosome), 'r') as f:
            output.append(list(map(float, f.readline().strip().split())))
    return np.array(output).sum(axis=0, keepdims=True).T


def get_tissue_parameter_snps(tissue: str, ancestry: str) -> List[str]:
    output = []
    for chromosome in range(1, 23):
        with open(tissue_m_path(tissue, ancestry, chromosome), 'r') as f:
            output.append(list(map(float, f.readline().strip().split())))
    return np.array(output).sum(axis=0, keepdims=True).T


def get_baseline_ld_score(ancestry: str) -> npt.NDArray:
    output = []
    for chromosome in range(1, 23):
        with gzip.open(baseline_path(ancestry, chromosome), 'rt') as f:
            _ = f.readline()  # unused header
            for line in f:
                output.append(list(map(float, line.strip().split('\t')[3:])))
    return np.array(output)


def get_tissue_ld_score(tissue: str, ancestry: str) -> npt.NDArray:
    output = []
    for chromosome in range(1, 23):
        with gzip.open(tissue_path(tissue, ancestry, chromosome), 'rt') as f:
            _ = f.readline()  # unused header
            for line in f:
                output.append(list(map(float, line.strip().split('\t')[3:])))
    return np.array(output)

def save_baseline_data(ancestry: str, ld: npt.NDArray, variables: List, parameter_snps: npt.NDArray) -> None:
    os.makedirs(f'inputs/{ancestry}/baseline/', exist_ok=True)
    np.save(f'inputs/{ancestry}/baseline/baseline_ld.{ancestry}.npy', ld)
    with open(f'inputs/{ancestry}/baseline/baseline_variables.{ancestry}.txt', 'w') as f:
        f.write('\t'.join(variables))
    np.save(f'inputs/{ancestry}/baseline/baseline_parameter_snps.{ancestry}.npy', parameter_snps)


def save_tissue_data(tissue: str, ancestry: str, ld: npt.NDArray, variables: List, parameter_snps: npt.NDArray) -> None:
    os.makedirs(f'inputs/{ancestry}/tissue/', exist_ok=True)
    np.save(f'inputs/{ancestry}/tissue/tissue_ld.{tissue}.{ancestry}.npy', ld)
    with open(f'inputs/{ancestry}/tissue/tissue_variables.{tissue}.{ancestry}.txt', 'w') as f:
        f.write('\t'.join(variables))
    np.save(f'inputs/{ancestry}/tissue/tissue_parameter_snps.{tissue}.{ancestry}.npy', parameter_snps)


def main():
    for ancestry in ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']:
        baseline_ld = get_baseline_ld_score(ancestry)
        baseline_variables = get_baseline_variables(ancestry)
        baseline_parameter_snps = get_baseline_parameter_snps(ancestry)
        save_baseline_data(ancestry, baseline_ld, baseline_variables, baseline_parameter_snps)

        for tissue_path in glob.glob(f'tissue/{ancestry}/*'):
            tissue = re.findall(f'tissue/{ancestry}/(.*)', tissue_path)[0]
            tissue_ld = get_tissue_ld_score(tissue, ancestry)
            tissue_variables = get_tissue_variables(tissue, ancestry)
            tissue_parameter_snps = get_tissue_parameter_snps(tissue, ancestry)
            save_tissue_data(tissue, ancestry, tissue_ld, tissue_variables, tissue_parameter_snps)


if __name__ == '__main__':
    main()
