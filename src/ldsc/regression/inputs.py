import numpy as np
from numpy import typing as npt
from typing import List


def baseline_ld_path(data_path: str, ancestry: str) -> str:
    return f'{data_path}/baseline/baseline_ld.{ancestry}.npy'


def baseline_variable_path(data_path: str, ancestry: str) -> str:
    return f'{data_path}/baseline/baseline_variables.{ancestry}.txt'


def baseline_parameter_snps_path(data_path: str, ancestry: str) -> str:
    return f'{data_path}/baseline/baseline_parameter_snps.{ancestry}.npy'


def tissue_ld_path(data_path: str, tissue: str, ancestry: str) -> str:
    return f'{data_path}/tissue/tissue_ld.{tissue}.{ancestry}.npy'


def tissue_variable_path(data_path: str, tissue: str, ancestry: str) -> str:
    return f'{data_path}/tissue/tissue_variables.{tissue}.{ancestry}.txt'


def tissue_parameter_snps_path(data_path: str, tissue: str, ancestry: str) -> str:
    return f'{data_path}/tissue/tissue_parameter_snps.{tissue}.{ancestry}.npy'


def overlap_path(data_path: str, ancestry: str) -> str:
    return f'{data_path}/overlap/overlap.baseline.{ancestry}.npy'


def overlap_tissue_path(data_path: str, tissue: str, ancestry: str) -> str:
    return f'{data_path}/overlap/overlap.{tissue}.{ancestry}.npy'


def overlap_tissue_baseline_path(data_path: str, tissue: str, ancestry: str) -> str:
    return f'{data_path}/overlap/overlap.baseline.{tissue}.{ancestry}.npy'


def get_baseline_ld(data_path: str, ancestry: str) -> npt.NDArray:
    return np.load(baseline_ld_path(data_path, ancestry))


def get_baseline_variables(data_path: str, ancestry: str) -> List[str]:
    with open(baseline_variable_path(data_path, ancestry), 'r') as f:
        return f.readline().strip().split('\t')


def get_baseline_parameter_snps(data_path: str, ancestry: str) -> npt.NDArray:
    return np.load(baseline_parameter_snps_path(data_path, ancestry))


def get_tissue_ld(data_path: str, tissue: str, ancestry: str) -> npt.NDArray:
    return np.load(tissue_ld_path(data_path, tissue, ancestry))


def get_tissue_variables(data_path: str, tissue: str, ancestry: str) -> List[str]:
    with open(tissue_variable_path(data_path, tissue, ancestry), 'r') as f:
        return f.readline().strip().split('\t')


def get_tissue_parameter_snps(data_path: str, tissue: str, ancestry: str) -> npt.NDArray:
    return np.load(tissue_parameter_snps_path(data_path, tissue, ancestry))


def get_overlap(data_path: str, tissue: str, ancestry: str) -> npt.NDArray:
    baseline_overlap = np.load(overlap_path(data_path, ancestry))
    tissue_overlap = np.load(overlap_tissue_path(data_path, tissue, ancestry))
    baseline_tissue_overlap = np.load(overlap_tissue_baseline_path(data_path, tissue, ancestry))
    return np.vstack((
        np.hstack((baseline_overlap, baseline_tissue_overlap)),
        np.hstack((baseline_tissue_overlap.T, tissue_overlap))
    ))
