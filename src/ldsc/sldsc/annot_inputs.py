import gzip
import numpy as np
from numpy import typing as npt
import re
from typing import Dict, List
from zipfile import ZipFile

import annotation


def get_input_zip_path(data_path: str, ancestry: str) -> str:
    return f'{data_path}/inputs/annotation_inputs.{ancestry}.zip'


def get_version(data_path: str, ancestry: str) -> str:
    with ZipFile(get_input_zip_path(data_path, ancestry)) as input_zip:
        with input_zip.open('version', 'r') as f:
            return f.readline().decode().strip()


def baseline_annot_path() -> str:
    return f'baseline/baseline_annot.npy'


def baseline_ld_path() -> str:
    return f'baseline/baseline_ld.npy'


def baseline_variable_path() -> str:
    return f'baseline/baseline_variables.txt'


def baseline_parameter_snps_path() -> str:
    return f'baseline/baseline_parameter_snps.npy'


def frq_path(data_path: str, ancestry: str, chromosome: int) -> str:
    return f'{data_path}/frq/{ancestry}/chr.{chromosome}.frq'


def weights_path(phenotype: str) -> str:
    return f'phenotypes/{phenotype}/weights.npy'


def sample_size_path(phenotype: str) -> str:
    return f'phenotypes/{phenotype}/sample_size.npy'


def y_path(phenotype: str) -> str:
    return f'phenotypes/{phenotype}/y.npy'


def idxs_path(phenotype: str) -> str:
    return f'phenotypes/{phenotype}/idxs.npy'


def phenotype_map_path(data_path: str) -> str:
    return f'{data_path}/phenotype_map/phenotype_map.tsv'


def get_all_phenotypes(data_path: str, ancestry: str) -> List[str]:
    with ZipFile(get_input_zip_path(data_path, ancestry)) as input_zip:
        return [phenotype for file in input_zip.filelist
                for phenotype in re.findall(r'phenotypes/([^/]*)/y.npy', file.filename)]


def get_weights(data_path: str, ancestry: str, phenotype: str) -> npt.NDArray:
    with ZipFile(get_input_zip_path(data_path, ancestry)) as input_zip:
        with input_zip.open(weights_path(phenotype), 'r') as f:
            return np.load(f)


def get_sample_size(data_path: str, ancestry: str, phenotype: str) -> npt.NDArray:
    with ZipFile(get_input_zip_path(data_path, ancestry)) as input_zip:
        with input_zip.open(sample_size_path(phenotype), 'r') as f:
            return np.load(f)


def get_y(data_path: str, ancestry: str, phenotype: str) -> npt.NDArray:
    with ZipFile(get_input_zip_path(data_path, ancestry)) as input_zip:
        with input_zip.open(y_path(phenotype), 'r') as f:
            return np.load(f)


def get_idxs(data_path: str, ancestry: str, phenotype: str) -> npt.NDArray:
    with ZipFile(get_input_zip_path(data_path, ancestry)) as input_zip:
        with input_zip.open(idxs_path(phenotype), 'r') as f:
            return np.load(f)


def get_baseline_ld(data_path: str, ancestry: str) -> npt.NDArray:
    with ZipFile(get_input_zip_path(data_path, ancestry)) as input_zip:
        with input_zip.open(baseline_ld_path(), 'r') as f:
            return np.load(f)


def get_baseline_variables(data_path: str, ancestry: str) -> List[str]:
    with ZipFile(get_input_zip_path(data_path, ancestry)) as input_zip:
        with input_zip.open(baseline_variable_path(), 'r') as f:
            return f.readline().decode().strip().split('\t')


def get_baseline_parameter_snps(data_path: str, ancestry: str) -> npt.NDArray:
    with ZipFile(get_input_zip_path(data_path, ancestry)) as input_zip:
        with input_zip.open(baseline_parameter_snps_path(), 'r') as f:
            return np.load(f)


def get_overlap(input_path: str, data_path: str, ancestry: str) -> npt.NDArray:
    with ZipFile(get_input_zip_path(input_path, ancestry)) as input_zip:
        with input_zip.open(baseline_annot_path(), 'r') as f:
            baseline_annot = np.load(f)
    out = []
    for chromosome in range(1, 23):
        with gzip.open(annotation.annotation_path(data_path, chromosome), 'rt') as annot_file, \
                open(frq_path(input_path, ancestry, chromosome), 'r') as frq_file:
            annot_file.readline(), frq_file.readline()
            for annot_line, frq_line in zip(annot_file, frq_file):
                frq = float(frq_line.strip().split('\t')[4])
                if 0.05 < frq < 0.95:
                    values = annot_line.strip().split('\t')
                    out.append(list(map(float, values)))
    np_array = np.hstack((baseline_annot, np.array(out)))
    return np_array.T.dot(np_array)

def get_phenotype_map(input_path: str) -> Dict[str, str]:
    out = {}
    with open(phenotype_map_path(input_path), 'r') as f:
        for line in f:
            phenotype, name = line.strip().split('\t')
            out[phenotype] = name
    return out
