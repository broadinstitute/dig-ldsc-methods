import gzip
import os
import subprocess
from typing import List, Dict, Optional

var_id_columns = ['chromosome', 'position']
input_path = os.environ.get('INPUT_PATH')
s3_path = os.environ.get('S3_BUCKET')


def check_snpmap(genome_build: str) -> None:
    if not os.path.exists(f'{input_path}/snpmap/sumstats.{genome_build}.snpmap'):
        cmd = f'./bootstrap/snpmap.bootstrap.sh {s3_path} {input_path} {genome_build}'
        subprocess.check_call(cmd, shell=True)


def get_rs_map(genome_build: str) -> Dict:
    check_snpmap(genome_build)
    rs_map = {}
    with open(f'{input_path}/snpmap/sumstats.{genome_build}.snpmap', 'r') as f:
        for line in f:
            chromosome, position, rs_id = line.strip().split('\t')
            rs_map[(chromosome, position)] = rs_id
    return rs_map


def get_p_value(line: Dict, col_map: Dict) -> float:
    return float(line[col_map['pValue']])


def get_n(line: Dict, col_map: Dict, effective_n: Optional[float]) -> float:
    if effective_n is not None:
        return effective_n
    else:
        return float(line[col_map['n']])


def valid_line(line: Dict, col_map: Dict, effective_n: Optional[float]) -> bool:
    return all([line.get(col_map[column]) is not None for column in var_id_columns]) and \
        col_map['pValue'] in line and line[col_map['pValue']] != '' and \
        (('n' in col_map and col_map['n'] in line) or effective_n is not None) and \
        0 < float(line[col_map['pValue']]) <= 1


def stream_to_data(file_path: str, rs_map: Dict, metadata: Dict) -> (List, Dict):
    out = []
    counts = {'all': 0, 'flipped': 0, 'error': 0}
    effective_n = metadata.get('effective_n')
    col_map = metadata['col_map']
    separator = metadata['separator']
    with gzip.open(file_path, 'rt') as f_in:
        header = f_in.readline().strip().split(separator)
        for json_string in f_in:
            line = dict(zip(header, json_string.strip().split(separator)))
            counts['all'] += 1
            try:
                if valid_line(line, col_map, effective_n):
                    chromosome, position = (line[col_map[c]] for c in var_id_columns)
                    if (chromosome, position) in rs_map:
                        rs_id = rs_map[(chromosome, position)]
                        p_value = float(line[col_map['pValue']])
                        n = get_n(line, col_map, effective_n)
                        out.append((rs_id, p_value, n))
            except ValueError:  # pValue, beta, or n not a value that can be converted to a float, skip
                counts['error'] += 1
    counts['final'] = len(out)
    return out, counts


def save_to_file(data_path: str, data: List) -> None:
    os.makedirs(f'{data_path}/magma/sumstats/', exist_ok=True)
    out_file = f'{data_path}/magma/sumstats/magma.sumstats.csv'
    with open(out_file, 'w') as f:
        f.write('SNP\tP\tN\n')
        for rs_id, p, n in data:
            f.write('{}\t{}\t{}\n'.format(rs_id, p, n))


def main(data_path: str, metadata: Dict) -> Dict:
    file = metadata['file']
    genome_build = metadata['genome_build']

    rs_map = get_rs_map(genome_build)
    data, counts = stream_to_data(f'{data_path}/raw/{file}', rs_map, metadata)
    metadata['counts'] = counts
    save_to_file(data_path, data)
    return metadata
