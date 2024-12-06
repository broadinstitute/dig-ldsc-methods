import numpy as np
import pytest
from src.ldsc.sldsc import xtx_xty


def test_weights() -> None:
    g1000_ld = np.array([[0.1, 0.4, 0.7], [0.2, 0.5, 0.8], [0.3, 0.6, 0.9]])
    sample_size = np.array([[1000], [2000], [3000]])
    chisq = np.array([[1.5], [2.0], [2.5]])
    parameter_snps = np.array([[500], [1000], [1500]])
    l_hm3 = np.array([[24 * 24 / 22 / 22], [36 * 36 / 31 / 31], [72 * 72 / 43 / 43]])
    weight = xtx_xty.get_weight(g1000_ld, l_hm3, sample_size, chisq, parameter_snps)
    assert weight[0][0] == pytest.approx(1 / 2)
    assert weight[1][0] == pytest.approx(1 / 3)
    assert weight[2][0] == pytest.approx(1 / 6)


def test_xtx() -> None:
    separators = [0, 2, 4]
    x1 = np.array([[0.1, 0.2, 0.3, 0.4], [0.5, 0.6, 0.7, 0.8], [0.9, 1.0, 1.1, 1.2], [1.3, 1.4, 1.5, 1.6]])
    x2 = np.array([[0.7, 1.1, 1.3, 1.7], [1.2, 0.8, 1.5, 1.9], [1.4, 1.6, 0.9, 2.1], [1.8, 2.0, 2.2, 1.0]])
    xtx_blocks = xtx_xty.get_general_xtx(x1, x2, separators)
    assert xtx_blocks.shape == (2, 4, 4)

    block_1_addition = [0.19, 0.19, 0.28, 0.36]
    block_1_start = [0.67, 0.51, 0.88, 1.12]
    for i in range(4):
        for j in range(len(xtx_blocks[0, i, :])):
            assert xtx_blocks[0, i, j] == pytest.approx(block_1_start[j] + block_1_addition[j] * i)

    block_2_addition = [0.32, 0.36, 0.31, 0.31]
    block_2_start = [3.6, 4.04, 3.67, 3.19]
    for i in range(4):
        for j in range(len(xtx_blocks[1, i, :])):
            assert xtx_blocks[1, i, j] == pytest.approx(block_2_start[j] + block_2_addition[j] * i)


def test_x() -> None:
    ld_matrix = np.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6], [0.7, 0.8, 0.9]])
    weights = np.array([[1.0], [0.9], [0.8]])
    n = np.array([[1000], [2000], [3000]])
    x = xtx_xty.get_x(ld_matrix, weights, n)
    assert x.shape == (3, 3)

    target = np.array([[0.05, 0.1, 0.15], [0.36, 0.45, 0.54], [0.84, 0.96, 1.08]])
    for i in range(len(target)):
        for j in range(len(target[i])):
            assert target[i][j] == pytest.approx(x[i][j])


def test_separators():
    snps = 10
    max_blocks = 2
    separators = xtx_xty.get_separators(snps, max_blocks)
    assert separators == [0, 5, 10]

    snps = 11
    max_blocks = 2
    separators = xtx_xty.get_separators(snps, max_blocks)
    assert separators == [0, 5, 11]

    snps = 10
    max_blocks = 3
    separators = xtx_xty.get_separators(snps, max_blocks)
    assert separators == [0, 3, 6, 10]

    snps = 10
    max_blocks = 11
    separators = xtx_xty.get_separators(snps, max_blocks)
    assert separators == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
