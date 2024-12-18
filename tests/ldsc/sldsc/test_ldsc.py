import numpy as np
from numpy import typing as npt
import pytest
from typing import List
from src.ldsc.sldsc import ldsc


def get_xtx_xty(xys: List[List]) -> (npt.NDArray, npt.NDArray):
    x = np.array([xy[:-1] + [1.0] for xy in xys])
    y = np.array([[xy[-1]] for xy in xys])
    return x.T.dot(x), x.T.dot(y)


def test_per_snp_heritability() -> None:
    xtx1, xty1 = get_xtx_xty([[1.0, 1.0, 5.0], [1.0, 2.0, 8.0], [2.0, 1.0, 7.0]])
    xtx2, xty2 = get_xtx_xty([[1.0, 1.0, 10.0], [1.0, 2.0, 16.0], [2.0, 1.0, 14.0]])
    assert ldsc.get_per_snp_heritability(np.array([xtx1]), np.array([xty1]), 1.0)[0][0] == pytest.approx(2.0)
    assert ldsc.get_per_snp_heritability(np.array([xtx1]), np.array([xty1]), 1.0)[1][0] == pytest.approx(3.0)
    assert ldsc.get_per_snp_heritability(np.array([xtx1]), np.array([xty1]), 2.0)[0][0] == pytest.approx(1.0)
    assert ldsc.get_per_snp_heritability(np.array([xtx1]), np.array([xty1]), 2.0)[1][0] == pytest.approx(1.5)
    assert ldsc.get_per_snp_heritability(np.array([xtx2]), np.array([xty2]), 2.0)[0][0] == pytest.approx(2.0)
    assert ldsc.get_per_snp_heritability(np.array([xtx2]), np.array([xty2]), 2.0)[1][0] == pytest.approx(3.0)
    assert ldsc.get_per_snp_heritability(np.array([xtx1, xtx2]), np.array([xty1, xty2]), 2.0)[0][0] == pytest.approx(1.5)
    assert ldsc.get_per_snp_heritability(np.array([xtx1, xtx2]), np.array([xty1, xty2]), 2.0)[1][0] == pytest.approx(2.25)


def test_biased_blocks_per_snp_heritability() -> None:
    xtx1, xty1 = get_xtx_xty([[1.0, 1.0, 5.0], [1.0, 2.0, 8.0], [2.0, 1.0, 7.0]])
    xtx2, xty2 = get_xtx_xty([[1.0, 1.0, 10.0], [1.0, 2.0, 16.0], [2.0, 1.0, 14.0]])
    for mean_sample_size in [1.0, 2.0]:
        blocks = ldsc.get_biased_block_per_snp_heritability(np.array([xtx1, xtx2]), np.array([xty1, xty2]), mean_sample_size)
        block1 = ldsc.get_per_snp_heritability(np.array([xtx1]), np.array([xty1]), mean_sample_size)
        block2 = ldsc.get_per_snp_heritability(np.array([xtx2]), np.array([xty2]), mean_sample_size)
        assert blocks.shape[0] == 2  # num blocks
        assert blocks.shape[1] == xtx1.shape[0] - 1  # remove intercept
        assert blocks.shape[2] == 1  # always 1
        assert all([blocks[0][i][j] == block2[i][j] for i in range(blocks.shape[1]) for j in range(blocks.shape[2])])  # combined - block1 = block2
        assert all([blocks[1][i][j] == block1[i][j] for i in range(blocks.shape[1]) for j in range(blocks.shape[2])])  # combined - block1 = block2


def test_unbiased_bloack_values() -> None:
    xtx1, xty1 = get_xtx_xty([[1.0, 1.0, 5.0], [1.0, 2.0, 8.0], [2.0, 1.0, 7.0]])
    xtx2, xty2 = get_xtx_xty([[1.0, 1.0, 10.0], [1.0, 2.0, 16.0], [2.0, 1.0, 14.0]])
    values = ldsc.get_per_snp_heritability(np.array([xtx1, xtx2]), np.array([xty1, xty2]), 2.0)
    biased_block_values = ldsc.get_biased_block_per_snp_heritability(np.array([xtx1, xtx2]), np.array([xty1, xty2]), 2.0)
    unbiased_block_values = ldsc.get_unbiased_block_values(values, biased_block_values)
    assert unbiased_block_values[0][0][0] == pytest.approx(1.0)
    assert unbiased_block_values[0][1][0] == pytest.approx(1.5)
    assert unbiased_block_values[1][0][0] == pytest.approx(2.0)
    assert unbiased_block_values[1][1][0] == pytest.approx(3.0)


def test_jackknife_covariance() -> None:
    xtx1, xty1 = get_xtx_xty([[1.0, 1.0, 5.0], [1.0, 2.0, 8.0], [2.0, 1.0, 7.0]])
    xtx2, xty2 = get_xtx_xty([[1.0, 1.0, 10.0], [1.0, 2.0, 16.0], [2.0, 1.0, 14.0]])
    values = ldsc.get_per_snp_heritability(np.array([xtx1, xtx2]), np.array([xty1, xty2]), 2.0)
    biased_block_values = ldsc.get_biased_block_per_snp_heritability(np.array([xtx1, xtx2]), np.array([xty1, xty2]), 2.0)
    unbiased_block_values = ldsc.get_unbiased_block_values(values, biased_block_values)
    cov = ldsc.jackknife_covariance(unbiased_block_values)
    assert cov[0][0] == pytest.approx(0.5 * 0.5)
    assert cov[0][1] == pytest.approx(0.5 * 0.75)
    assert cov[1][0] == pytest.approx(0.75 * 0.5)
    assert cov[1][1] == pytest.approx(0.75 * 0.75)


def test_heritability_proportion() -> None:
    xtx1, xty1 = get_xtx_xty([[1.0, 1.0, 5.0], [1.0, 2.0, 8.0], [2.0, 1.0, 7.0]])
    xtx2, xty2 = get_xtx_xty([[1.0, 1.0, 10.0], [1.0, 2.0, 16.0], [2.0, 1.0, 14.0]])
    per_snp_heritability = ldsc.get_per_snp_heritability(np.array([xtx1, xtx2]), np.array([xty1, xty2]), 2.0)
    heritability_proportion = ldsc.get_heritability_proportion(per_snp_heritability, np.array([[1.0], [1.0]]))
    assert heritability_proportion[0][0] == pytest.approx(0.4)
    assert heritability_proportion[1][0] == pytest.approx(0.6)

    heritability_proportion = ldsc.get_heritability_proportion(per_snp_heritability, np.array([[2.0], [4.0]]))
    assert heritability_proportion[0][0] == pytest.approx(0.25)
    assert heritability_proportion[1][0] == pytest.approx(0.75)


def test_biased_block_hertiability_proportion() -> None:
    xtx1, xty1 = get_xtx_xty([[1.0, 1.0, 5.0], [1.0, 2.0, 8.0], [2.0, 1.0, 7.0]])
    xtx2, xty2 = get_xtx_xty([[1.0, 1.0, 10.0], [1.0, 2.0, 16.0], [2.0, 1.0, 14.0]])
    biased_block_values = ldsc.get_biased_block_per_snp_heritability(np.array([xtx1, xtx2]), np.array([xty1, xty2]), 2.0)
    biased_block_heritability_proportion = ldsc.get_biased_block_heritability_proportion(biased_block_values, np.array([[2.0], [4.0]]))
    assert biased_block_heritability_proportion[0][0][0] == pytest.approx(0.25)
    assert biased_block_heritability_proportion[0][1][0] == pytest.approx(0.75)
    assert biased_block_heritability_proportion[1][0][0] == pytest.approx(0.25)
    assert biased_block_heritability_proportion[1][1][0] == pytest.approx(0.75)


def test_corrected_se() -> None:
    xtx1, xty1 = get_xtx_xty([[1.0, 1.0, 5.0], [1.0, 2.0, 8.0], [2.0, 1.0, 7.0]])
    xtx2, xty2 = get_xtx_xty([[1.0, 1.0, 10.0], [1.0, 2.0, 16.0], [2.0, 1.0, 14.0]])
    correction_matrix = np.array([[1.0, 2.0], [3.0, 4.0]])

    values = ldsc.get_per_snp_heritability(np.array([xtx1, xtx2]), np.array([xty1, xty2]), 2.0)
    biased_block_values = ldsc.get_biased_block_per_snp_heritability(np.array([xtx1, xtx2]), np.array([xty1, xty2]), 2.0)
    unbiased_block_values = ldsc.get_unbiased_block_values(values, biased_block_values)
    cov = ldsc.jackknife_covariance(unbiased_block_values)
    corrected_cov = ldsc.get_corrected_se(cov, correction_matrix)

    a = 0.5 * 0.5
    b = 0.5 * 0.75
    c = 0.75 * 0.75
    assert corrected_cov[0][0] == pytest.approx(np.sqrt(1 * a + 2 * b + 2 * b + 4 * c))
    assert corrected_cov[1][0] == pytest.approx(np.sqrt(9 * a + 12 * b + 12 * b + 16 * c))


def test_integrated_h2() -> None:
    xtx1, xty1 = get_xtx_xty([[1.0, 1.0, 5.0], [1.0, 2.0, 8.0], [2.0, 1.0, 7.0]])
    xtx2, xty2 = get_xtx_xty([[1.0, 1.0, 10.0], [1.0, 2.0, 16.0], [2.0, 1.0, 14.0]])
    xtx, xty = np.array([xtx1, xtx2]), np.array([xty1, xty2])
    overlap_matrix = np.array([[100, 50], [50, 100]])
    parameter_snps = np.array([[2.0], [4.0]])
    total_snps = 0.32
    mean_sample_size = 2.0
    h2 = ldsc.get_h2(xtx, xty, overlap_matrix, parameter_snps, total_snps, mean_sample_size)

    assert pytest.approx(h2[0]['expHeritability']) == 6.25
    assert pytest.approx(h2[0]['heritability']) == 21.875
    assert pytest.approx(h2[0]['enrichment']) == 3.5
    assert pytest.approx(h2[0]['pValue']) == 0.095466

    assert pytest.approx(h2[1]['expHeritability']) == 12.5
    assert pytest.approx(h2[1]['heritability']) == 25.0
    assert pytest.approx(h2[1]['enrichment']) == 2.0
    assert pytest.approx(h2[1]['pValue']) == 0.095466
