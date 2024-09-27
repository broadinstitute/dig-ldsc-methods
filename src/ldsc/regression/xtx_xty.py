import numpy as np
from numpy import typing as npt
from typing import List


def get_weight(g1000_ld: npt.NDArray, l_hm3: npt.NDArray, sample_size: npt.NDArray, chisq: npt.NDArray, parameter_snps: npt.NDArray) -> npt.NDArray:
    g1000_ld_sum = np.sum(g1000_ld, axis=1, keepdims=True)
    tau = (np.mean(chisq) - 1) / np.mean(np.multiply(g1000_ld_sum, sample_size))
    print(tau)

    safe_g1000_ld_sum = np.fmax(g1000_ld_sum, 1.0)
    safe_tau = min(max(tau, 0.0), 1.0 / np.sum(parameter_snps))
    safe_l_hm3 = np.fmax(l_hm3, 1.0)

    heteroskedasticity_weight = 1.0 / (1 + safe_tau * sample_size * safe_g1000_ld_sum)**2
    over_counting_weight = 1.0 / safe_l_hm3
    print(np.mean(over_counting_weight), np.mean(heteroskedasticity_weight))
    unnormalized_weight = np.sqrt(heteroskedasticity_weight * over_counting_weight)
    return unnormalized_weight / np.sum(unnormalized_weight)


def get_general_xtx(x1: npt.NDArray, x2: npt.NDArray, separators: List[int]) -> npt.NDArray:
    blocks = len(separators) - 1
    xtx_blocks = np.zeros((blocks, x1.shape[1], x2.shape[1]))
    for i in range(blocks):
        xtx_blocks[i, :, :] = x1[separators[i]:separators[i + 1], :].T.dot(x2[separators[i]:separators[i + 1], :])
    return xtx_blocks


def get_x(ld_matrix: npt.NDArray, weights: npt.NDArray, n: npt.NDArray) -> npt.NDArray:
    conditioned_ld_matrix = np.multiply(n, ld_matrix) / np.mean(n)  # keep condition number low
    return np.multiply(conditioned_ld_matrix, weights)


def get_intercept(size: int, weights: npt.NDArray) -> npt.NDArray:
    intercept_col = np.ones((size, 1))
    return np.multiply(intercept_col, weights)


def get_y(chisq: npt.NDArray, weights: npt.NDArray) -> npt.NDArray:
    return np.multiply(chisq, weights)


def get_xtx(x: npt.NDArray, separators: List[int]) -> npt.NDArray:
    return get_general_xtx(x, x, separators)


def get_xty(x: npt.NDArray, y: npt.NDArray, separators: List[int]) -> npt.NDArray:
    return get_general_xtx(x, y, separators)


def get_separators(snps: int, max_blocks: int) -> List[int]:
    return list(map(int, np.floor(np.linspace(0, snps, min(max_blocks, snps) + 1))))
