import time
import gzip
import os
import subprocess
from typing import Dict, Iterator, List, Tuple

input_path = os.environ.get('INPUT_PATH')
s3_path = os.environ.get('S3_BUCKET')

def check_envvars() -> None:
    assert input_path is not None
    assert s3_path is not None


def g1000_path(ancestry: str, chromosome: str) -> str:
    return f'{input_path}/g1000/{ancestry}/chr{chromosome}.bim'


def check_g1000(ancestry: str) -> None:
    if not os.path.exists(g1000_path(ancestry, '1')):
        cmd = f'./bootstrap/g1000.bootstrap.sh {s3_path} {input_path} {ancestry}'
        subprocess.check_call(cmd, shell=True)


def annotation_path(data_path: str) -> str:
    return f'{data_path}/raw/annotation.tsv'


def get_annotation_data(data_path: str, file_chromosome: str) -> List[Tuple[int, int]]:
    data = []
    with open(annotation_path(data_path), 'r') as f:
        for line in f:
            chromosome, start, end = line.strip().split('\t', 2)
            if '\t' in end:
                end, _ = end.split('\t', 1)
            if chromosome == file_chromosome:
                data.append((int(start), int(end)))
    return data


def get_g1000_data(ancestry: str, chromosome: str) -> List[int]:
    data = []
    with open(g1000_path(ancestry, chromosome), 'r') as f:
        for line in f:
            _, _, _, position, _ = line.strip().split('\t', 4)
            data.append(int(position))
    return data


def next_range(iter_range: Iterator[Tuple[int, int]], max_val: int):
    return next(iter_range, (max_val, max_val))


def write_annot(data_path: str, chromosome: str, range_data: List[Tuple[int, int]], g1000_data: List) -> None:
    os.makedirs(f'{data_path}/sldsc/annot-ld/', exist_ok=True)
    out_file = f'{data_path}/sldsc/annot-ld/ld.{chromosome}.annot.gz'
    with gzip.open(out_file, 'wt') as f_out:
        f_out.write('ANNOT\n')
        iter_range = iter(range_data)
        curr_start, curr_end = next_range(iter_range, g1000_data[-1])
        for curr_position in iter(g1000_data):
            while curr_end < curr_position:
                curr_start, curr_end = next_range(iter_range, g1000_data[-1])
            if curr_start < curr_position <= curr_end:
                f_out.write('1\n')
            else:
                f_out.write('0\n')


def run_chromosome(data_path: str, ancestry: str, chromosome: str) -> None:
    range_data = get_annotation_data(data_path, chromosome)
    g1000_data = get_g1000_data(ancestry, chromosome)
    write_annot(data_path, chromosome, range_data, g1000_data)


def annotation(data_path: str, metadata: Dict) -> None:
    ancestry = metadata['ancestry']

    check_g1000(ancestry)
    tot_time = 0
    for chromosome in range(1, 23):
        t = time.time()
        run_chromosome(data_path, ancestry, str(chromosome))
        print(f'Chromosome Run Time ({chromosome}): {time.time() - t}')
        tot_time += time.time() - t
    print(f'Full Annotation Run Time: {tot_time}')
