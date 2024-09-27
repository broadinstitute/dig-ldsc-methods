import gzip
import numpy as np
from numpy import typing as npt


def dataset_path(data_path: str, dataset: str) -> str:
    return f'{data_path}/dataset/{dataset}.sumstats.gz'


def load_sumstats(data_path: str, dataset: str) -> (npt.NDArray, npt.NDArray, npt.NDArray):
    idxs = []
    z = []
    n = []
    snps = []
    with gzip.open(dataset_path(data_path, dataset), 'rt') as f:
        _ = f.readline()  # header
        for i, line in enumerate(f):
            split_line = line.strip().split('\t')
            if len(split_line) == 3:
                snps.append(split_line[0])
                z.append(float(split_line[1]))
                n.append(float(split_line[2]))
                idxs.append(i)
    return np.array([z]).T**2, np.array([n]).T, np.array(idxs)


def get_filter_idxs(chisq: npt.NDArray, n: npt.NDArray) -> npt.NDArray:
    chisq_max = max(0.001 * np.max(n), 80)
    return chisq[:, 0] < chisq_max


def filter_sumstats(chisq: npt.NDArray, n: npt.NDArray, idxs: npt.NDArray) -> (npt.NDArray, npt.NDArray, npt.NDArray):
    filter_idxs = get_filter_idxs(chisq, n)
    return chisq[filter_idxs, :], n[filter_idxs], idxs[filter_idxs]
