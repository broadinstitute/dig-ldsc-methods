import gzip
import numpy as np
from numpy import typing as npt
from typing import List


def annotation_path(data_path: str, chromosome: int) -> str:
    return f'{data_path}/sldsc/annot-ld/ld.{chromosome}.annot.gz'


def annotation_ld_path(data_path: str, chromosome: int) -> str:
    return f'{data_path}/sldsc/annot-ld/ld.{chromosome}.l2.ldscore.gz'


def annotation_m_path(data_path: str, chromosome: int) -> str:
    return f'{data_path}/sldsc/annot-ld/ld.{chromosome}.l2.M_5_50'


def get_ld(data_path: str) -> npt.NDArray:
    output = []
    for chromosome in range(1, 23):
        with gzip.open(annotation_ld_path(data_path, chromosome), 'rt') as f:
            _ = f.readline()
            for line in f:
                output.append(float(line.strip().split('\t')[3]))
    return np.array([output]).T


def get_parameter_snps(ancestry: str) -> List[str]:
    output = []
    for chromosome in range(1, 23):
        with open(annotation_m_path(ancestry, chromosome), 'r') as f:
            output.append(list(map(float, f.readline().strip().split())))
    return np.array(output).sum(axis=0, keepdims=True).T
