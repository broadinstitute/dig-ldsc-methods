import gzip
import numpy as np
from numpy import typing as npt


def weights_path(data_path: str, ancestry: str, chromosome: int) -> str:
    return f'{data_path}/weights/{ancestry}/weights.{chromosome}.l2.ldscore.gz'


def get_input_weights(data_path: str, ancestry: str) -> npt.NDArray:
    output = []
    for chromosome in range(1, 23):
        with gzip.open(weights_path(data_path, ancestry, chromosome), 'rt') as f:
            _ = f.readline()
            line = f.readline().strip().split('\t')
            while len(line) > 1:
                output.append(float(line[3]))
                line = f.readline().strip().split('\t')
    return np.array([output]).T
