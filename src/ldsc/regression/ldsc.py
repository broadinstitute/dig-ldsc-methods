import numpy as np
from numpy import typing as npt
from scipy.stats import t as tdist
from typing import Dict, List, Optional


# In:  xtx (blocks, parameters + intercept, parameters + intercept)
#      xty (blocks, parameters + intercept, 1)
#      mean_sample_size (scalar)
# Out: per_snp_heritability (parameters, 1) NOTE: intercept removed
def get_per_snp_heritability(xtx: npt.NDArray, xty: npt.NDArray, mean_sample_size: float) -> npt.NDArray:
    return np.linalg.solve(np.sum(xtx, axis=0), np.sum(xty, axis=0))[:-1, :] / mean_sample_size


# In:  xtx (blocks, parameters + intercept, parameters + intercept)
#      xty (blocks, parameters + intercept, 1)
#      mean_sample_size (scalar)
# Out: biased_block_per_snp_heritability (blocks, parameters, 1) NOTE: intercept removed
def get_biased_block_per_snp_heritability(xtx: npt.NDArray, xty: npt.NDArray, mean_sample_size: float) -> npt.NDArray:
    blocks, parameters, _ = xtx.shape
    biased_block_per_snp_heritability = np.zeros((blocks, parameters, 1))
    xty_tot = np.sum(xty, axis=0)
    xtx_tot = np.sum(xtx, axis=0)
    for j in range(blocks):
        # undo conditioning
        biased_block_per_snp_heritability[j, :, :] = \
            np.linalg.solve(xtx_tot - xtx[j], xty_tot - xty[j]) / mean_sample_size
    return biased_block_per_snp_heritability[:, :-1, :]  # remove intercept


# In:  values (parameters, 1)
#      biased_block_values (blocks, parameters, 1)
# Out: unbiased_block_values (blocks, parameters, 1)
def get_unbiased_block_values(values: npt.NDArray, biased_block_values: npt.NDArray) -> npt.NDArray:
    blocks = biased_block_values.shape[0]
    return blocks * values - (blocks - 1) * biased_block_values


# In:  unbiased_block_values (blocks, parameters, 1)
# Out: covariance (parameters, parameters)
def jackknife_covariance(unbiased_block_values: npt.NDArray) -> npt.NDArray:
    return np.cov(unbiased_block_values[:, :, 0].T, ddof=1) / unbiased_block_values.shape[0]


# In:  values (parameters, 1)
#      biased_block_values (blocks, parameters, 1)
# Out: covariance (parameters, parameters)
def get_covariance(values: npt.NDArray, biased_block_values: npt.NDArray) -> npt.NDArray:
    unbiased_block_values = get_unbiased_block_values(values, biased_block_values)
    return jackknife_covariance(unbiased_block_values)


# In:  per_snp_heritability (parameters, 1)
#      parameter_snps (parameters, 1)
# Out: heritability_proportion (parameters, 1)
def get_heritability_proportion(per_snp_heritability: npt.NDArray, parameter_snps: npt.NDArray) -> npt.NDArray:
    heritability = parameter_snps * per_snp_heritability
    return heritability / np.sum(heritability)


# In:  biased_block_per_snp_heritability (blocks, parameters, 1)
#      parameter_snps (parameters, 1)
# Out: biased_block_heritability_proportion (blocks, parameters, 1)
def get_biased_block_heritability_proportion(biased_block_per_snp_heritability: npt.NDArray, parameter_snps: npt.NDArray) -> npt.NDArray:
    biased_block_heritability_proportion = np.zeros(biased_block_per_snp_heritability.shape)
    for j in range(biased_block_per_snp_heritability.shape[0]):
        mult = np.multiply(parameter_snps, biased_block_per_snp_heritability[j, :, :])
        biased_block_heritability_proportion[j, :, :] = mult / np.sum(mult)
    return biased_block_heritability_proportion


# In:  covariance_matrix (parameters, parameters)
#      correction_matrix (parameters, parameters)
# Out: corrected_se (parameters, 1)
def get_corrected_se(covariance_matrix: npt.NDArray, correction_matrix: npt.NDArray) -> npt.NDArray:
    variances = correction_matrix.dot(covariance_matrix).dot(correction_matrix.T).diagonal()
    return np.vstack(np.sqrt(np.maximum(0, variances)))


# Student-T survival function (1 - cdf) with blocks = dof
def get_p_value(mu: float, se: float, blocks: int) -> Optional[float]:
    if se != 0:
        return 2 * tdist.sf(abs(mu / se), blocks)


def get_h2(xtx: npt.NDArray,  # (blocks, parameters + intercept, parameters + intercept)
           xty: npt.NDArray,  # (blocks, parameters + intercept, 1)
           overlap_matrix: npt.NDArray,  # (parameters, parameters)
           parameter_snps: npt.NDArray,  # (parameters, 1)
           total_snps: float,  # (scalar)
           mean_sample_size: float  # (scalar)
           ) -> List[Dict]:  # list(json)
    blocks = xtx.shape[0]

    per_snp_heritability = get_per_snp_heritability(xtx, xty, mean_sample_size)
    biased_block_per_snp_heritability = get_biased_block_per_snp_heritability(xtx, xty, mean_sample_size)
    per_snp_heritability_covariance = get_covariance(per_snp_heritability, biased_block_per_snp_heritability)

    heritability_proportion = get_heritability_proportion(per_snp_heritability, parameter_snps)
    biased_block_heritability_proportion = get_biased_block_heritability_proportion(biased_block_per_snp_heritability, parameter_snps)
    heritability_proportion_covariance = get_covariance(heritability_proportion, biased_block_heritability_proportion)

    # adjust based on overlap (equation 3: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4626285/)
    attenuation_correction = np.multiply(overlap_matrix, 1 / parameter_snps.T)
    corrected_heritability_proportion = attenuation_correction.dot(heritability_proportion)
    corrected_heritability_proportion_se = get_corrected_se(heritability_proportion_covariance, attenuation_correction)

    expected_heritability_proportion = parameter_snps / total_snps
    enrichment = corrected_heritability_proportion / expected_heritability_proportion
    enrichment_se = corrected_heritability_proportion_se / expected_heritability_proportion

    parameter_snp_difference = total_snps - parameter_snps
    attenuation_correction = np.multiply(overlap_matrix, total_snps / parameter_snps / parameter_snp_difference) - \
                             parameter_snps.T / parameter_snp_difference
    corrected_per_snp_heritability = attenuation_correction.dot(per_snp_heritability)
    corrected_per_snp_heritability_se = get_corrected_se(per_snp_heritability_covariance, attenuation_correction)

    p_value = [
        get_p_value(corrected_per_snp_heritability[i][0], corrected_per_snp_heritability_se[i][0], blocks)
        for i in range(len(corrected_per_snp_heritability_se))
    ]

    return [{
        'expHeritability': expected_heritability_proportion[i][0],
        'heritability': corrected_heritability_proportion[i][0],
        'heritabilitySE': corrected_heritability_proportion_se[i][0],
        'enrichment': enrichment[i][0],
        'enrichmentSE': enrichment_se[i][0],
        'pValue': p_value[i] if p_value[i] is not None else 'NA'
    } for i in range(len(enrichment))]
