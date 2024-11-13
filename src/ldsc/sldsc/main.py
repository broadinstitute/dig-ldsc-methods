import argparse
import json
import numpy as np
import os
import subprocess
from typing import Dict

import inputs, ldsc, sumstats, weights, xtx_xty

MAX_BLOCKS = 200
input_path = os.environ['INPUT_PATH']


def check_sldsc_inputs(ancestry: str) -> None:
    if not os.path.exists(inputs.get_input_zip_path(input_path, ancestry)):
        subprocess.check_call(f'./bootstrap/sldsc.bootstrap.sh {input_path} {ancestry}', shell=True)


def check_weights(ancestry: str) -> None:
    if not os.path.exists(weights.weights_path(input_path, ancestry, 1)):
        subprocess.check_call(f'./bootstrap/weights.bootstrap.sh {input_path} {ancestry}', shell=True)


def get_metadata(data_path: str) -> Dict:
    with open(f'{data_path}/sumstats/metadata', 'r') as f:
        metadata = json.load(f)
    return metadata


def save_data(output: Dict, metadata: Dict, data_path: str) -> None:
    os.makedirs(f'{data_path}/sldsc/', exist_ok=True)
    for variable_type, data in output.items():
        file = f'{data_path}/sldsc/{variable_type}.output.tsv'
        with open(file, 'w') as f:
            for line in data:
                f.write(line)
    with open(f'{data_path}/sldsc/metadata', 'w') as f:
        metadata['version'] = inputs.get_version(input_path, metadata['ancestry'])
        json.dump(metadata, f)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir', default=None, required=True, type=str)
    data_path = parser.parse_args().dir

    metadata = get_metadata(data_path)
    ancestry = metadata['ancestry']

    check_sldsc_inputs(ancestry)
    check_weights(ancestry)

    chisq, sample_size, idxs = sumstats.load_sumstats(data_path)
    chisq, sample_size, idxs = sumstats.filter_sumstats(chisq, sample_size, idxs)
    mean_sample_size = float(np.mean(sample_size))

    baseline_ld = inputs.get_baseline_ld(input_path, ancestry)
    baseline_variables = inputs.get_baseline_variables(input_path, ancestry)
    baseline_parameter_snps = inputs.get_baseline_parameter_snps(input_path, ancestry)
    input_weights = weights.get_input_weights(input_path, ancestry)

    sumstats_weights = input_weights[idxs]
    baseline_sumstats_ld = baseline_ld[idxs, :]
    separators = xtx_xty.get_separators(baseline_sumstats_ld.shape[0], MAX_BLOCKS)

    baseline_weights = xtx_xty.get_weight(baseline_sumstats_ld, sumstats_weights, sample_size, chisq, baseline_parameter_snps)
    baseline_x = xtx_xty.get_x(baseline_sumstats_ld, baseline_weights, sample_size)
    intercept = xtx_xty.get_intercept(baseline_sumstats_ld.shape[0], baseline_weights)
    y = xtx_xty.get_y(chisq, baseline_weights)

    output = {}
    for tissue in inputs.get_all_tissues(input_path, ancestry):
        overlap_matrix = inputs.get_overlap(input_path, tissue, ancestry)
        total_snps = overlap_matrix[0][0]

        variables = baseline_variables + inputs.get_tissue_variables(input_path, tissue, ancestry)
        tissue_sumstats_ld = inputs.get_tissue_ld(input_path, tissue, ancestry)[idxs, :]
        parameter_snps = np.vstack((baseline_parameter_snps, inputs.get_tissue_parameter_snps(input_path, tissue, ancestry)))

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
    save_data(output, metadata, data_path)


if __name__ == '__main__':
    main()
