import argparse
import gzip
import json
import math
import numpy as np
import os
from scipy.stats import chi2
from typing import List, Dict, Optional

var_id_columns = ['chromosome', 'position', 'reference', 'alt']
input_path = os.environ['INPUT_PATH']


def get_metadata(data_path: str) -> Dict:
    with open(f'{data_path}/raw/metadata', 'r') as f:
        metadata = json.load(f)
    return metadata


def get_var_to_rs_map(ancestry: str, genome_build: str, build_type: str) -> Dict:
    var_to_rs = {}
    with open(f'{input_path}/snpmap/sumstats.{build_type}.{genome_build}.{ancestry}.snpmap', 'r') as f:
        for full_row in f.readlines():
            var_id, rs_id = full_row.strip().split('\t')
            var_to_rs[var_id] = rs_id
    return var_to_rs


def get_p_value(line: Dict, col_map: Dict) -> float:
    return float(line[col_map['pValue']])


def get_beta(line: Dict, col_map: Dict) -> float:
    if 'beta' in col_map:
        return float(line[col_map['beta']])
    else:
        return math.log(float(line[col_map['oddsRatio']]))


def get_n(line: Dict, col_map: Dict, effective_n: Optional[float]) -> float:
    if effective_n is not None:
        return effective_n
    else:
        return float(line[col_map['n']])


def p_to_z(p: float, beta: float) -> float:
    return np.sqrt(chi2.isf(p, 1)) * (-1)**(beta < 0)


def valid_line(line: Dict, col_map: Dict, effective_n: Optional[float]) -> bool:
    return all([line.get(col_map[column]) is not None for column in var_id_columns]) and \
           col_map['pValue'] in line and \
           (('beta' in col_map and col_map['beta'] in line) or ('oddsRatio' in col_map and col_map['oddsRatio'] in line)) and \
           (('n' in col_map and col_map['n'] in line) or effective_n is not None) and \
           0 < float(line[col_map['pValue']]) <= 1


def stream_to_data(data_path: str, file: str, ancestry: str, genome_build: str, effective_n: Optional[float], col_map: Dict, separator: str) -> List:
    var_to_rs_map = get_var_to_rs_map(ancestry, genome_build, 'standard')
    var_to_rs_flipped = get_var_to_rs_map(ancestry, genome_build, 'flipped')
    out = []
    count = 0
    error_count = 0
    flip_count = 0
    with gzip.open(f'{data_path}/raw/{file}', 'rt') as f_in:
        header = f_in.readline().strip().split(separator)
        for json_string in f_in:
            line = dict(zip(header, json_string.strip().split(separator)))
            count += 1
            if valid_line(line, col_map, effective_n):
                chromosome, position, reference, alt = (line[col_map[c]] for c in var_id_columns)
                var_id = f'{chromosome}:{position}:{reference.upper()}:{alt.upper()}'
                if var_id in var_to_rs_map:
                    try:
                        p_value = float(line[col_map['pValue']])
                        beta = get_beta(line, col_map)
                        n = get_n(line, col_map, effective_n)
                        out.append((var_to_rs_map[var_id], p_to_z(p_value, beta), n))
                    except ValueError:  # pValue, beta, or n not a value that can be converted to a float, skip
                        error_count += 1
                elif var_id in var_to_rs_flipped:
                    try:
                        p_value = float(line[col_map['pValue']])
                        beta = get_beta(line, col_map) * -1  # flip
                        n = get_n(line, col_map, effective_n)
                        out.append((var_to_rs_flipped[var_id], p_to_z(p_value, beta), n))
                        flip_count += 1
                    except ValueError:  # pValue, beta, or n not a value that can be converted to a float, skip
                        error_count += 1
    print('{} SNPs translated from a total of {} ({} skipped, {} flipped)'.format(len(out), count, error_count, flip_count))
    return out


def filter_data_to_dict(data: List) -> Dict:
    N90 = np.quantile([a[2] for a in data], 0.9)
    return {d[0]: (d[1], d[2]) for d in data if d[2] >= N90 / 1.5}


def save_to_file(data_path: str, ancestry: str, data: Dict, metadata: Dict) -> None:
    os.makedirs(f'{data_path}/sumstats/', exist_ok=True)
    out_file = f'{data_path}/sumstats/sldsc.sumstats.gz'
    with gzip.open(out_file, 'wt') as f:
        f.write('SNP\tZ\tN\n')
        for rs_id in ld_rs_iter(ancestry):
            if rs_id in data:
                f.write('{}\t{}\t{}\n'.format(rs_id, round(data[rs_id][0], 3), data[rs_id][1]))
            else:
                f.write(f'{rs_id}\t\t\n')
    with open(f'{data_path}/sumstats/metadata', 'w') as f:
        json.dump(metadata, f)


def ld_rs_iter(ancestry: str) -> str:
    for CHR in range(1, 23):
        with gzip.open(f'{input_path}/weights/{ancestry}/weights.{CHR}.l2.ldscore.gz', 'rt') as f:
            _ = f.readline()
            for line in f:
                yield line.strip().split('\t')[1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir', default=None, required=True, type=str)
    data_path = parser.parse_args().dir

    metadata = get_metadata(data_path)
    file = metadata['file']
    ancestry = metadata['ancestry']
    effective_n = metadata.get('effective_n')
    col_map = metadata['col_map']
    separator = metadata['separator']
    genome_build = metadata['genome_build']

    data = stream_to_data(data_path, file, ancestry, genome_build, effective_n, col_map, separator)
    if len(data) > 0:
        data_dict = filter_data_to_dict(data)
        save_to_file(data_path, ancestry, data_dict, metadata)


if __name__ == '__main__':
    main()
