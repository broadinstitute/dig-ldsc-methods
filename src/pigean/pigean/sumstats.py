import gzip
import os
from typing import Dict, Optional

var_id_columns = ['chromosome', 'position']
input_path = os.environ.get('INPUT_PATH')
s3_path = os.environ.get('S3_BUCKET')


def open_text(path: str):
    """Open an uploaded GWAS file for reading whether or not it is gzipped.

    The upload UI accepts both and the stored metadata does not say which
    was sent, so sniff the gzip magic bytes rather than trusting the name.
    """
    with open(path, 'rb') as f:
        is_gzip = f.read(2) == b'\x1f\x8b'
    return gzip.open(path, 'rt') if is_gzip else open(path, 'r')


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


def stream_to_sumstats(data_path: str, file: str, metadata: Dict) -> Dict:
    os.makedirs(f'{data_path}/pigean/sumstats/', exist_ok=True)
    counts = {'all': 0, 'translated': 0, 'skipped': 0, 'error': 0}
    effective_n = metadata.get('effective_n')
    col_map = metadata['col_map']
    separator = metadata['separator']
    with gzip.open(f'{data_path}/pigean/sumstats/pigean.sumstats.gz', 'wt') as f_out:
        f_out.write('CHROM\tPOS\tP\tN\n')
        with open_text(f'{data_path}/raw/{file}') as f_in:
            header = f_in.readline().strip().split(separator)
            for json_string in f_in:
                line = dict(zip(header, json_string.strip().split(separator)))
                counts['all'] += 1
                try:
                    if valid_line(line, col_map, effective_n):
                        chromosome, position = (line[col_map[c]] for c in var_id_columns)
                        p_value = float(line[col_map['pValue']])
                        n = get_n(line, col_map, effective_n)
                        f_out.write('{}\t{}\t{}\t{}\n'.format(chromosome, position, p_value, n))
                        counts['translated'] += 1
                    else:
                        counts['skipped'] += 1
                except ValueError:  # pValue, beta, or n not a value that can be converted to a float, skip
                    counts['error'] += 1
    return counts


def main(data_path: str, metadata: Dict) -> Dict:
    file = metadata['file']
    counts = stream_to_sumstats(data_path, file, metadata)
    metadata['counts'] = counts
    return metadata
