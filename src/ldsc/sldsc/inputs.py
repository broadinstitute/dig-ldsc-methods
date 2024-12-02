import numpy as np
from numpy import typing as npt
import re
from typing import List
from zipfile import ZipFile


def get_input_zip_path(data_path: str, ancestry: str) -> str:
    return f'{data_path}/inputs/sldsc_inputs.{ancestry}.zip'


def get_version(data_path: str, ancestry: str) -> str:
    with ZipFile(get_input_zip_path(data_path, ancestry)) as input_zip:
        with input_zip.open('version', 'r') as f:
            return f.readline().decode().strip()


def get_all_tissues(data_path: str, ancestry: str) -> List[str]:
    tissues = set()
    with ZipFile(get_input_zip_path(data_path, ancestry)) as input_zip:
        for tissue_path in input_zip.namelist():
            maybe_tissue = re.findall(f'tissue/tissue_ld.(.*).{ancestry}.npy', tissue_path)
            if len(maybe_tissue) > 0:
                tissues |= {maybe_tissue[0]}
    return sorted(tissues)


def baseline_ld_path(ancestry: str) -> str:
    return f'baseline/baseline_ld.{ancestry}.npy'


def baseline_variable_path(ancestry: str) -> str:
    return f'baseline/baseline_variables.{ancestry}.txt'


def baseline_parameter_snps_path(ancestry: str) -> str:
    return f'baseline/baseline_parameter_snps.{ancestry}.npy'


def tissue_ld_path(tissue: str, ancestry: str) -> str:
    return f'tissue/tissue_ld.{tissue}.{ancestry}.npy'


def tissue_variable_path(tissue: str, ancestry: str) -> str:
    return f'tissue/tissue_variables.{tissue}.{ancestry}.txt'


def tissue_parameter_snps_path(tissue: str, ancestry: str) -> str:
    return f'tissue/tissue_parameter_snps.{tissue}.{ancestry}.npy'


def overlap_path(ancestry: str) -> str:
    return f'overlap/overlap.baseline.{ancestry}.npy'


def overlap_tissue_path(tissue: str, ancestry: str) -> str:
    return f'overlap/overlap.{tissue}.{ancestry}.npy'


def overlap_tissue_baseline_path(tissue: str, ancestry: str) -> str:
    return f'overlap/overlap.baseline.{tissue}.{ancestry}.npy'


def get_baseline_ld(data_path: str, ancestry: str) -> npt.NDArray:
    with ZipFile(get_input_zip_path(data_path, ancestry)) as input_zip:
        with input_zip.open(baseline_ld_path(ancestry), 'r') as f:
            return np.load(f)


def get_baseline_variables(data_path: str, ancestry: str) -> List[str]:
    with ZipFile(get_input_zip_path(data_path, ancestry)) as input_zip:
        with input_zip.open(baseline_variable_path(ancestry), 'r') as f:
            return f.readline().decode().strip().split('\t')


def get_baseline_parameter_snps(data_path: str, ancestry: str) -> npt.NDArray:
    with ZipFile(get_input_zip_path(data_path, ancestry)) as input_zip:
        with input_zip.open(baseline_parameter_snps_path(ancestry), 'r') as f:
            return np.load(f)


def get_tissue_ld(data_path: str, tissue: str, ancestry: str) -> npt.NDArray:
    with ZipFile(get_input_zip_path(data_path, ancestry)) as input_zip:
        with input_zip.open(tissue_ld_path(tissue, ancestry), 'r') as f:
            return np.load(f)


def get_tissue_variables(data_path: str, tissue: str, ancestry: str) -> List[str]:
    with ZipFile(get_input_zip_path(data_path, ancestry)) as input_zip:
        with input_zip.open(tissue_variable_path(tissue, ancestry), 'r') as f:
            return f.readline().decode().strip().split('\t')


def get_tissue_parameter_snps(data_path: str, tissue: str, ancestry: str) -> npt.NDArray:
    with ZipFile(get_input_zip_path(data_path, ancestry)) as input_zip:
        with input_zip.open(tissue_parameter_snps_path(tissue, ancestry), 'r') as f:
            return np.load(f)


def get_overlap(data_path: str, tissue: str, ancestry: str) -> npt.NDArray:
    with ZipFile(get_input_zip_path(data_path, ancestry)) as input_zip:
        with input_zip.open(overlap_path(ancestry), 'r') as f:
            baseline_overlap = np.load(f)
        with input_zip.open(overlap_tissue_path(tissue, ancestry), 'r') as f:
            tissue_overlap = np.load(f)
        with input_zip.open(overlap_tissue_baseline_path(tissue, ancestry), 'r') as f:
            baseline_tissue_overlap = np.load(f)
        return np.vstack((
            np.hstack((baseline_overlap, baseline_tissue_overlap)),
            np.hstack((baseline_tissue_overlap.T, tissue_overlap))
        ))
