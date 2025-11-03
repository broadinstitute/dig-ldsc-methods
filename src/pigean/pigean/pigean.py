import gzip
import json
import os
from typing import Callable, Dict, List
import subprocess

input_path = os.environ.get('INPUT_PATH')
s3_path = os.environ.get('S3_BUCKET')


def get_pigean_models_path(input_path: str) -> str:
    return f'{input_path}/models/aws_pigean_models_s3.json'


def check_pigean_models() -> None:
    if not os.path.exists(get_pigean_models_path(input_path)):
        subprocess.check_call(f'./bootstrap/pigean.bootstrap.sh {s3_path} {input_path}', shell=True)


def get_model_data() -> (Dict, Dict):
    with open(f'{input_path}/models/aws_pigean_models_s3.json', 'r') as f:
        models = json.load(f)
    return ({model['name']: model for model in models['models']},
            {gene_set['name']: gene_set for gene_set in models['gene_sets']})


def file_name(data_path: str, trait_type: str) -> str:
    if trait_type == 'sumstats':
        return f'{data_path}/pigean/sumstats/pigean.sumstats.gz'
    elif trait_type == 'gene_lists':
        return f'{data_path}/pigean/sumstats/gene_list.tsv'
    elif trait_type == 'exomes':
        return f'{data_path}/pigean/sumstats/exomes.sumstats.gz'
    else:
        raise ValueError(f'Invalid trait_type: {trait_type}')


def get_gene_sets(model: str) -> List[str]:
    models, gene_sets = get_model_data()
    model_info = models[model]
    inputs = []
    p_infs = []
    for gene_set in model_info['gene_sets']:
        gene_set_info = gene_sets[gene_set]
        if gene_set_info['type'] == 'set':
            inputs += ['--X-in', f'{input_path}/gene_sets/{gene_set_info["file"]}']
        else:
            inputs += ['--X-list', f'{input_path}/gene_lists/{gene_set_info["name"]}/{gene_set_info["file"]}']
        p_infs += ['--p-noninf', str(gene_set_info['p-inf'])]
    if len(inputs) > 0:
        return inputs + p_infs
    else:
        raise Exception(f'Invalid model {model}')


def input_type_command(data_path: str, input_type: str) -> List[str]:
    if input_type == 'sumstats':
        return ['--gwas-in', file_name(data_path, input_type),
                '--gwas-chrom-col', 'CHROM',
                '--gwas-pos-col', 'POS',
                '--gwas-p-col', 'P',
                '--gwas-n-col', 'N'
                ]
    elif input_type == 'gene_lists':
        return [
            '--positive-controls-in', file_name(data_path, input_type),
            '--positive-controls-id-col', '1',
            '--positive-controls-prob-col', '2',
            '--positive-controls-no-header', 'True',
            '--positive-controls-all-in', f'{input_path}/misc/refGene_hg19_TSS.subset.loc',
            '--positive-controls-all-no-header', 'True',
            '--positive-controls-all-id-col', '1'
        ]
    elif input_type == 'exomes':
        return [
            '--exomes-in', file_name(data_path, input_type),
            '--exomes-gene-col', 'Gene',
            '--exomes-p-col', 'P-value',
            '--exomes-beta-col', 'Effect'
        ]


def base_cmd(data_path: str) -> List[str]:
    return [
        'python3', f'{input_path}/methods/priors.py', 'gibbs',
        '--first-for-hyper',
        '--sigma-power', '-2',
        '--gwas-detect-high-power', '100',
        '--gwas-detect-low-power', '10',
        '--num-chains', '10',
        '--num-chains-betas', '4',
        '--max-num-iter', '500',
        '--filter-gene-set-p', '0.005',
        '--max-num-gene-sets', '5000',
        '--gene-filter-value', '1.0',
        '--gene-set-filter-value', '0.01',
        '--s2g-normalize-values', '0.95',
        '--update-hyper', 'none',
        '--gene-loc-file', f'{input_path}/misc/NCBI37.3.plink.gene.loc',
        '--gene-map-in', f'{input_path}/misc/portal_gencode.gene.map',
        '--gene-loc-file-huge', f'{input_path}/misc/refGene_hg19_TSS.subset.loc',
        '--exons-loc-file-huge', f'{input_path}/misc/NCBI37.3.plink.gene.exons.loc',
        '--gene-stats-out', f'{data_path}/pigean/pigean/gene_stats.out',
        '--gene-set-stats-out', f'{data_path}/pigean/pigean/gene_set_stats.out',
        '--gene-gene-set-stats-out', f'{data_path}/pigean/pigean/gene_gene_set_stats.out',
        '--gene-effectors-out', f'{data_path}/pigean/pigean/gene_effector.out'
    ]


def run_pigean(data_path: str, input_type: str, model: str) -> None:
    os.makedirs(f'{data_path}/pigean/pigean/', exist_ok=True)
    cmd = base_cmd(data_path) + input_type_command(data_path, input_type) + get_gene_sets(model)
    subprocess.check_call(cmd)


def make_option(value: str) -> str:
    return value if value != 'NA' else 'null'


def translate_gene_stats(json_line: Dict, model: str) -> str:
    combined = make_option(json_line["combined"])
    huge_score = json_line["huge_score_gwas"] if 'huge_score_gwas' in json_line else json_line["positive_control"]
    if combined != 'null':
        return f'{{"gene": "{json_line["Gene"]}", ' \
               f'"prior": {make_option(json_line["prior"])}, ' \
               f'"combined": {combined}, ' \
               f'"huge_score": {make_option(huge_score)}, ' \
               f'"log_bf": {make_option(json_line["log_bf"])}, ' \
               f'"n": {make_option(json_line["N"])}, ' \
               f'"model": "{model}"}}\n'


def translate_gene_set_stats(json_line: Dict, model: str):
    beta = make_option(json_line["beta"])
    beta_uncorrected = make_option(json_line["beta_uncorrected"])
    if beta != 'null' and beta_uncorrected != 'null' and float(beta_uncorrected) != 0.0:
        return f'{{"gene_set": "{json_line["Gene_Set"]}", ' \
               f'"source": "{json_line["label"]}", ' \
               f'"beta": {beta}, ' \
               f'"beta_uncorrected": {beta_uncorrected}, ' \
               f'"n": {make_option(json_line["N"])}, ' \
               f'"model": "{model}"}}\n'


def get_gene_set_stats_maps(data_path: str) -> (Dict, Dict):
    beta_uncorrected_map = {}
    source_map = {}
    with open(f'{data_path}/pigean/pigean/gene_set_stats.out', 'r') as f_in:
        header = f_in.readline().strip().split('\t')
        for line in f_in:
            json_line = dict(zip(header, line.strip().split('\t')))
            source_map[json_line['Gene_Set']] = json_line['label']
            beta_uncorrected = make_option(json_line['beta_uncorrected'])
            if beta_uncorrected != 'null':
                beta_uncorrected_map[json_line['Gene_Set']] = json_line['beta_uncorrected']
    return beta_uncorrected_map, source_map


def get_translate_gene_gene_set_stats(data_path: str) -> Callable[[Dict, str], str]:
    beta_uncorrected_map, source_map = get_gene_set_stats_maps(data_path)
    def translate_gene_gene_set_stats(json_line, model):
        beta = make_option(json_line["beta"])
        combined = make_option(json_line["combined"])
        if beta is not None and combined is not None:
            beta_uncorrected = beta_uncorrected_map.get(json_line['gene_set'], '0.0')
            source = source_map[json_line['gene_set']]
            return f'{{"gene": "{json_line["Gene"]}", ' \
                   f'"gene_set": "{json_line["gene_set"]}", ' \
                   f'"source": "{source}", ' \
                   f'"prior": {make_option(json_line["prior"])}, ' \
                   f'"combined": {combined}, ' \
                   f'"beta": {beta}, ' \
                   f'"beta_uncorrected": {beta_uncorrected}, ' \
                   f'"log_bf": {make_option(json_line["log_bf"])}, ' \
                   f'"model": "{model}"}}\n'
    return translate_gene_gene_set_stats

def translate_gene_effector(json_line: Dict, model: str) -> str:
    return f'{{"lead_locus": "{json_line["Lead_locus"]}", ' \
           f'"p": "{json_line["P"]}", ' \
           f'"gene": {json_line["Gene"]}, ' \
           f'"cond_prob_total": {json_line["cond_prob_total"]}, ' \
           f'"cond_prob_signal": {json_line["cond_prob_signal"]}, ' \
           f'"cond_prob_prior": {json_line["cond_prob_prior"]}, ' \
           f'"cond_prob_huge": {json_line["cond_prob_huge"]}, ' \
           f'"combined_D": {json_line["combined_D"]}, ' \
           f'"model": "{model}"}}\n'


def convert(data_path: str, output_type: str, translate_fnc: Callable[[Dict, str], str], model) -> None:
    if os.path.exists(f'{data_path}/pigean/pigean/{output_type}.out'):
        with gzip.open(f'{data_path}/pigean/pigean/{output_type}.json.gz', 'wt') as f_out:
            with open(f'{data_path}/pigean/pigean/{output_type}.out', 'r') as f:
                header = f.readline().strip().split('\t')
                for line in f:
                    maybe_output = translate_fnc(dict(zip(header, line.strip().split('\t'))), model)
                    if maybe_output is not None:
                        f_out.write(maybe_output)


def main(data_path: str, metadata: Dict) -> Dict:
    check_pigean_models()
    run_pigean(data_path, metadata['input_type'], metadata['model'])
    convert(data_path, 'gene_stats', translate_gene_stats, metadata['model'])
    convert(data_path, 'gene_set_stats', translate_gene_set_stats, metadata['model'])
    convert(data_path, 'gene_gene_set_stats', get_translate_gene_gene_set_stats(data_path), metadata['model'])
    convert(data_path, 'gene_effector', translate_gene_effector, metadata['model'])
    return metadata
