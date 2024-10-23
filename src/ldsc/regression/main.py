import glob
import json
import numpy as np
import os
import re
from typing import Dict, List

import inputs, ldsc, sumstats, weights, xtx_xty

MAX_BLOCKS = 200
data_path = os.environ.get('DATA_PATH', '.')


def get_all_tissues(ancestry: str) -> List[str]:
    tissues = set()
    for tissue_path in glob.glob(f'{data_path}/inputs/{ancestry}/tissue/tissue_ld.*___*.{ancestry}.npy'):
        tissue = re.findall(f'{data_path}/inputs/{ancestry}/tissue/tissue_ld.(.*).{ancestry}.npy', tissue_path)[0]
        tissues |= {tissue}
    return sorted(tissues)


def get_metadata():
    with open('dataset/metadata', 'r') as f:
        metadata = json.load(f)
    return metadata


def save_data(output: Dict) -> None:
    os.makedirs(f'regression/', exist_ok=True)
    for variable_type, data in output.items():
        file = f'regression/{variable_type}.output.tsv'
        with open(file, 'w') as f:
            for line in data:
                f.write(line)


def main():
    metadata = get_metadata()
    ancestry = metadata['ancestry']

    chisq, sample_size, idxs = sumstats.load_sumstats(data_path)
    chisq, sample_size, idxs = sumstats.filter_sumstats(chisq, sample_size, idxs)
    mean_sample_size = float(np.mean(sample_size))

    baseline_ld = inputs.get_baseline_ld(data_path, ancestry)
    baseline_variables = inputs.get_baseline_variables(data_path, ancestry)
    baseline_parameter_snps = inputs.get_baseline_parameter_snps(data_path, ancestry)
    input_weights = weights.get_input_weights(data_path, ancestry)

    sumstats_weights = input_weights[idxs]
    baseline_sumstats_ld = baseline_ld[idxs, :]
    separators = xtx_xty.get_separators(baseline_sumstats_ld.shape[0], MAX_BLOCKS)

    baseline_weights = xtx_xty.get_weight(baseline_sumstats_ld, sumstats_weights, sample_size, chisq, baseline_parameter_snps)
    baseline_x = xtx_xty.get_x(baseline_sumstats_ld, baseline_weights, sample_size)
    intercept = xtx_xty.get_intercept(baseline_sumstats_ld.shape[0], baseline_weights)
    y = xtx_xty.get_y(chisq, baseline_weights)

    output = {}
    for tissue in get_all_tissues(ancestry):
        overlap_matrix = inputs.get_overlap(data_path, tissue, ancestry)
        total_snps = overlap_matrix[0][0]

        variables = baseline_variables + inputs.get_tissue_variables(data_path, tissue, ancestry)
        tissue_sumstats_ld = inputs.get_tissue_ld(data_path, tissue, ancestry)[idxs, :]
        parameter_snps = np.vstack((baseline_parameter_snps, inputs.get_tissue_parameter_snps(data_path, tissue, ancestry)))

        tissue_x = xtx_xty.get_x(tissue_sumstats_ld, baseline_weights, sample_size)

        x = np.hstack((baseline_x, tissue_x, intercept))
        xtx = xtx_xty.get_xtx(x, separators)
        xty = xtx_xty.get_xty(x, y, separators)

        values = ldsc.get_h2(xtx, xty, overlap_matrix, parameter_snps, total_snps, mean_sample_size)
        annotation, tissue_name = tissue.split('___')
        for variable, value in zip(variables, values):
            variable_type, variable_name = variable.split('___')
            if variable_type not in output:
                output[variable_type] = []
            line = f'{annotation}\t{tissue_name}\t{variable_name}\t{value["enrichment"]}\t{value["pValue"]}\n'
            output[variable_type].append(line)
    save_data(output)


if __name__ == '__main__':
    main()
