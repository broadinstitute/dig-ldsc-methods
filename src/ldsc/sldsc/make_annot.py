import time
import gzip
import os
import subprocess
from typing import Dict, Iterator, List, Tuple

input_path = os.environ.get('INPUT_PATH')
s3_path = os.environ.get('S3_BUCKET')

GENE_WINDOW = 50000

def check_envvars() -> None:
    assert input_path is not None
    assert s3_path is not None


def g1000_path(ancestry: str, chromosome: str) -> str:
    return f'{input_path}/g1000/{ancestry}/chr{chromosome}.bim'


def check_g1000(ancestry: str) -> None:
    if not os.path.exists(g1000_path(ancestry, '1')):
        cmd = f'./bootstrap/g1000.bootstrap.sh {s3_path} {input_path} {ancestry}'
        subprocess.check_call(cmd, shell=True)


def gene_loc_path() -> str:
    return f'{input_path}/gene_loc/gene.loc'


def check_gene_loc() -> None:
    if not os.path.exists(gene_loc_path()):
        cmd = f'./bootstrap/gene_loc.bootstrap.sh {s3_path} {input_path}'
        subprocess.check_call(cmd, shell=True)


def annotation_path(data_path: str, file: str) -> str:
    return f'{data_path}/raw/{file}'


def get_gene_loc_map():
    gene_loc_map = {}
    with open(gene_loc_path(), 'r') as f:
        for line in f:
            gene, chromosome, start, end = line.strip().split('\t')
            gene_loc_map[gene] = (chromosome, int(start), int(end))
    return gene_loc_map


def convert_gene_list_to_bed_file(data_path: str, file: str) -> str:
    gene_loc_map = get_gene_loc_map()
    with open(annotation_path(data_path, 'gene_list.bed'), 'w') as f_out:
        with open(annotation_path(data_path, file), 'r') as f:
            for line in f:
                gene = line.strip()
                if gene in gene_loc_map:
                    chromosome, start, end = gene_loc_map[gene]
                    f_out.write('{}\t{}\t{}\n'.format(
                        chromosome,
                        max(start - GENE_WINDOW, 1),
                        end + GENE_WINDOW
                    ))
    return 'gene_list.bed'


def get_annotation_data(data_path: str, file: str, file_chromosome: str) -> List[Tuple[int, int]]:
    data = []
    with open(annotation_path(data_path, file), 'r') as f:
        for line in f:
            split_line = line.strip().split('\t')
            if len(split_line) >= 3:
                chromosome, start, end = split_line[:3]
                if 'chr' in chromosome:
                    chromosome = chromosome.replace('chr', '')
                if chromosome == file_chromosome:
                    try:
                        data.append((int(start), int(end)))
                    except:
                        pass
    return sorted(data)


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


def run_chromosome(data_path: str, file: str, ancestry: str, chromosome: str) -> None:
    range_data = get_annotation_data(data_path, file, chromosome)
    g1000_data = get_g1000_data(ancestry, chromosome)
    write_annot(data_path, chromosome, range_data, g1000_data)


def annotation(data_path: str, metadata: Dict) -> None:
    ancestry = metadata['ancestry']
    file = metadata['file']
    file_type = metadata['file_type']

    if file_type == 'gene_list':
        check_gene_loc()
        file = convert_gene_list_to_bed_file(data_path, file)

    check_g1000(ancestry)
    tot_time = 0
    for chromosome in range(1, 23):
        t = time.time()
        run_chromosome(data_path, file, ancestry, str(chromosome))
        print(f'Chromosome Run Time ({chromosome}): {time.time() - t}')
        tot_time += time.time() - t
    print(f'Full Annotation Run Time: {tot_time}')
