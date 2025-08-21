import json
import numpy as np
import os
import subprocess
from typing import Dict

import annotation, annot_inputs, ldsc, xtx_xty

MAX_BLOCKS = 200
input_path = os.environ.get('INPUT_PATH')
s3_path = os.environ.get('S3_BUCKET')


def check_envvars():
    assert input_path is not None
    assert s3_path is not None


def check_inputs(ancestry: str) -> None:
    if not os.path.exists(annot_inputs.get_input_zip_path(input_path, ancestry)):
        subprocess.check_call(f'./bootstrap/annotation.bootstrap.sh {s3_path} {input_path} {ancestry}', shell=True)


def check_frq(ancestry: str) -> None:
    if not os.path.exists(annot_inputs.frq_path(input_path, ancestry, 1)):
        cmd = f'./bootstrap/frq.bootstrap.sh {s3_path} {input_path} {ancestry}'
        subprocess.check_call(cmd, shell=True)


def get_metadata(data_path: str) -> Dict:
    with open(f'{data_path}/sldsc/annot-ld/metadata', 'r') as f:
        metadata = json.load(f)
    return metadata


def save_data(output: Dict, metadata: Dict, data_path: str) -> None:
    os.makedirs(f'{data_path}/sldsc/annot-sldsc/', exist_ok=True)
    for variable_type, data in output.items():
        file = f'{data_path}/sldsc/annot-sldsc/{variable_type}.output.tsv'
        with open(file, 'w') as f:
            for line in data:
                f.write(line)
    with open(f'{data_path}/sldsc/annot-sldsc/metadata', 'w') as f:
        metadata['version'] = annot_inputs.get_version(input_path, metadata['ancestry'])
        json.dump(metadata, f)


def annot_sldsc(data_path):
    check_envvars()

    metadata = get_metadata(data_path)
    ancestry = metadata['ancestry']

    check_inputs(ancestry)
    check_frq(ancestry)

    overlap_matrix = annot_inputs.get_overlap(input_path, data_path, ancestry)
    total_snps = overlap_matrix[0][0]

    baseline_variables = annot_inputs.get_baseline_variables(input_path, ancestry)
    baseline_ld = annot_inputs.get_baseline_ld(input_path, ancestry)
    baseline_parameter_snps = annot_inputs.get_baseline_parameter_snps(input_path, ancestry)
    annotation_variable = ['custom___annotation']
    annotation_ld = annotation.get_ld(data_path)
    annotation_parameter_snps = annotation.get_parameter_snps(data_path)

    variables = baseline_variables + annotation_variable
    parameter_snps = np.vstack((baseline_parameter_snps, annotation_parameter_snps))

    output = {}
    phenotypes = annot_inputs.get_all_phenotypes(input_path, ancestry)
    for i, phenotype in enumerate(phenotypes):
        baseline_weights = annot_inputs.get_weights(input_path, ancestry, phenotype)
        sample_size = annot_inputs.get_sample_size(input_path, ancestry, phenotype)
        y = annot_inputs.get_y(input_path, ancestry, phenotype)
        idxs = annot_inputs.get_idxs(input_path, ancestry, phenotype)
        mean_sample_size = float(np.mean(sample_size))

        baseline_x = xtx_xty.get_x(baseline_ld[idxs, :], baseline_weights, sample_size)
        annotation_x = xtx_xty.get_x(annotation_ld[idxs, :], baseline_weights, sample_size)
        x = np.hstack((baseline_x, annotation_x, baseline_weights))

        separators = xtx_xty.get_separators(x.shape[0], MAX_BLOCKS)
        xtx = xtx_xty.get_xtx(x, separators)
        xty = xtx_xty.get_xty(x, y, separators)

        values = ldsc.get_h2(xtx, xty, overlap_matrix, parameter_snps, total_snps, mean_sample_size)
        for variable, value in zip(variables, values):
            variable_type, variable_name = variable.split('___')
            if variable_type not in output:
                output[variable_type] = []
            line = f'{phenotype}\t{ancestry}\t{variable_name}\t{value["enrichment"]}\t{value["pValue"]}\n'
            output[variable_type].append(line)
            if variable_type == 'custom':
                print(i, line.strip())
    save_data(output, metadata, data_path)
