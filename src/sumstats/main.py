import argparse
import gzip
import json
import numpy as np
from scipy.stats import chi2
import subprocess

data_path = '.'

def download(username, dataset):
    path = f's3://dig-ldsc-server/upload/{username}/{dataset}/'
    subprocess.check_call(f'aws s3 cp {path} dataset/ --recursive', shell=True)
    metadata = json.load('dataset/metadata')
    return metadata['file'], metadata['ancestry'], metadata['col_map']


def get_var_to_rs_map(ancestry):
    var_to_rs = {}
    with open(f'{data_path}/snplist/ldsc.{ancestry}.snplist', 'r') as f:
        for full_row in f.readlines():
            var_id, rs_id = full_row.strip().split('\t')
            var_to_rs[var_id] = rs_id
    return var_to_rs


def p_to_z(p, beta):
    return np.sqrt(chi2.isf(p, 1)) * (-1)**(beta < 0)


required_columns = ['chromosome', 'position', 'reference', 'alt', 'beta', 'pValue', 'n']
def valid_line(line, col_map):
    return all([line.get(col_map[c]) is not None for c in required_columns]) and 0 < line[col_map['pValue']] <= 1


var_id_columns = ['chromosome', 'position', 'reference', 'alt']
zn_columns = ['pValue', 'beta', 'n']
def stream_to_data(file, ancestry, col_map):
    var_to_rs_map = get_var_to_rs_map(ancestry)
    out = []
    with open(file, 'r') as f_in:
        for json_string in f_in:
            line = json.loads(json_string)
            if valid_line(line, col_map):
                chromosome, position, reference, alt = (line[col_map[c]] for c in var_id_columns)
                var_id = f'{chromosome}:{position}:{reference}:{alt}'
                if var_id in var_to_rs_map:
                    pValue, beta, n = (line[col_map[c]] for c in zn_columns)
                    out.append(var_to_rs_map[var_id], p_to_z(float(pValue), float(beta)), float(n))
    return out


def ld_rs_iter(ancestry):
    for CHR in range(1, 23):
        with gzip.open(f'{data_path}/weights/{ancestry}/weights.{CHR}.l2.ldscore.gz', 'rt') as f:
            _ = f.readline()
            for line in f:
                yield line.strip().split('\t')[1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--username', default=None, required=True, type=str)
    parser.add_argument('--dataset', default=None, required=True, type=str)
    args = parser.parse_args()

    file, ancestry, col_map = download(args.username, args.dataset)
    data = stream_to_data(file, ancestry, col_map)

    N90 = np.quantile([a[2] for a in data], 0.9)
    data = {d[0]: (d[1], d[2]) for d in data if d[2] >= N90 / 1.5}

    with gzip.open(f'data/{args.dataset}/{args.dataset}.sumstats.gz', 'wt') as f:
        f.write('SNP\tZ\tN\n')
        for rs_id in ld_rs_iter(ancestry):
            if rs_id in data:
                f.write('{}\t{}\t{}\n'.format(rs_id, *data[rs_id]))
            else:
                f.write('\t\t\n')


if __name__ == '__main__':
    main()
