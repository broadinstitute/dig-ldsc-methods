import argparse
import bitarray as ba
import gzip
import numpy as np
from numpy import typing as npt
import os
import shutil
import subprocess
from typing import List, Set

s3_bucket = 's3://dig-ldsc-server'
data_path = '.'

bedcode = {
    2: ba.bitarray('11'),
    1: ba.bitarray('01'),
    0: ba.bitarray('00')
}


def get_hapmap_set(chromosome: int) -> Set:
    with open(f'{data_path}/hapmap/hm.{chromosome}.snp', 'r') as f:
        return {line.strip() for line in f}


def get_cm(ancestry: str, chromosome: int) -> List:
    cm = []
    with open(f'g1000/{ancestry}/chr{chromosome}.bim', 'r') as f:
        for line in f:
            split_line = line.strip().split('\t')
            cm.append(float(split_line[2]))
    return cm


def get_hapmap_idxs(ancestry: str, chromosome: int, hapmap_set: Set) -> List:
    hapmap_idxs = []
    with open(f'g1000/{ancestry}/chr{chromosome}.bim', 'r') as f:
        for idx, line in enumerate(f):
            split_line = line.strip().split('\t')
            if split_line[1] in hapmap_set:
                hapmap_idxs.append(idx)
    return hapmap_idxs


def get_rs_positions(ancestry: str, chromosome: int, hapmap_set: Set) -> List:
    rs_positions = []
    with open(f'g1000/{ancestry}/chr{chromosome}.bim', 'r') as f:
        for idx, line in enumerate(f):
            split_line = line.strip().split('\t')
            if split_line[1] in hapmap_set:
                rs_positions.append((split_line[1], split_line[3]))
    return rs_positions


def download(username: str, dataset: str, ancestry: str) -> None:
    path = f'{s3_bucket}/userdata/{username}/annotation/{dataset}/ldsc/annotation'
    file = f'{dataset}.{ancestry}.zip'
    subprocess.check_call(f'aws s3 cp {path}/{file} data/annotation/', shell=True)
    subprocess.check_call(f'unzip data/annotation/{file} -d data/annotation/', shell=True)
    os.remove(f'data/annotation/{file}')


def get_annotation_idxs(dataset: str, ancestry: str, chromosome: int) -> List:
    with gzip.open(f'data/annotation/{dataset}.{ancestry}.{chromosome}.annot.gz', 'rt') as f:
        f.readline()  # header
        return [idx for idx, line in enumerate(f) if line == '1\n']


def check_bits(f, byte_count: int, target_string: str, error_str: str) -> None:
    file_bytes = ba.bitarray(endian="little")
    file_bytes.fromfile(f, byte_count)
    if file_bytes != ba.bitarray(target_string):
        raise IOError(error_str)


def get_genotype(ancestry: str, chromosome: str) -> ba.bitarray:
    with open(f'g1000/{ancestry}/chr{chromosome}.bed', 'rb') as f:
        check_bits(f, 2, '0011011011011000', 'Magic number from Plink .bed file not recognized')
        check_bits(f, 1, '10000000', 'Plink .bed file must be in default SNP-major mode')
        genotype = ba.bitarray(endian="little")
        genotype.fromfile(f)
        return genotype


def get_people(ancestry: str, chromosome: int) -> int:
    with open(f'g1000/{ancestry}/chr{chromosome}.fam', 'r') as f:
        return len(f.readlines())


def get_people_bits(people: int) -> int:
    extra_bits = (4 - people % 4) if people % 4 != 0 else 0
    return (people + extra_bits) * 2


def get_matrix(idxs: List, genotype: ba.bitarray, people_bits: int, people: int) -> npt.NDArray:
    raw_array = [genotype[people_bits * idx: people_bits * (idx + 1)].decode(bedcode) for idx in idxs]
    return np.array(raw_array)[:, :people]


def filter_heterogenous(matrix: npt.NDArray, idxs: List, people: int) -> (npt.NDArray, List):
    filter = np.sum(matrix == 1, 1) != people
    filtered_idxs = [idx for filter_in, idx in zip(filter, idxs) if filter_in]
    return matrix[filter, :], filtered_idxs


def normalize_matrix(matrix: npt.NDArray) -> npt.NDArray:
    mean = np.mean(matrix, 1, keepdims=True)
    std = np.std(matrix, 1, keepdims=True)
    ones_array = np.ones((1, matrix.shape[1]))
    return (matrix - mean.dot(ones_array)) / std.dot(ones_array)


def get_left_right_pairs(left_cm: List, right_cm: List, window_cm: float = 1.0) -> List:
    left_right_pairs = []
    idx = left_idx = right_idx = 0
    while idx < len(left_cm):
        while right_idx < len(right_cm) and right_cm[right_idx] - left_cm[idx] <= window_cm:
            right_idx += 1
        while left_idx < len(right_cm) and left_cm[idx] - right_cm[left_idx] > window_cm:
            left_idx += 1
        left_right_pairs.append((left_idx, right_idx - 1))
        idx += 1
    return left_right_pairs


def write_ldscore(dataset: str, ancestry: str, chromosome: int, ldscores: List, rs_positions: List) -> None:
    os.makedirs(f'data/ldscore/', exist_ok=True)
    with gzip.open(f'data/ldscore/{dataset}.{ancestry}.{chromosome}.l2.ldscore.gz', 'w') as f:
        f.write(b'CHR\tSNP\tBP\tL2\n')
        for (rs_id, position), ldscore in zip(rs_positions, ldscores):
            f.write(f'{chromosome}\t{rs_id}\t{position}\t{round(ldscore, 3)}\n'.encode())


def write_m(dataset: str, ancestry: str, chromosome: int, m: int, m_5_50: int) -> None:
    with open(f'data/ldscore/{dataset}.{ancestry}.{chromosome}.l2.M', 'w') as f:
        f.write(f'{m}\n')
    with open(f'data/ldscore/{dataset}.{ancestry}.{chromosome}.l2.M_5_50', 'w') as f:
        f.write(f'{m_5_50}\n')


def run(dataset: str, ancestry: str, chromosome: str) -> None:
    hapmap_set = get_hapmap_set(chromosome)
    genotype = get_genotype(ancestry, chromosome)
    people = get_people(ancestry, chromosome)
    people_bits = get_people_bits(people)

    hapmap_idxs = get_hapmap_idxs(ancestry, chromosome, hapmap_set)
    annotation_idxs = get_annotation_idxs(dataset, ancestry, chromosome)

    hapmap_matrix = get_matrix(hapmap_idxs, genotype, people_bits, people)
    annotation_matrix = get_matrix(annotation_idxs, genotype, people_bits, people)
    annotation_matrix, annotation_idxs = filter_heterogenous(annotation_matrix, annotation_idxs, people)

    hapmap_matrix = normalize_matrix(hapmap_matrix)
    annotation_matrix = normalize_matrix(annotation_matrix)

    cm = get_cm(ancestry, chromosome)
    hapmap_cm = np.array(cm)[hapmap_idxs]
    annotation_cm = np.array(cm)[annotation_idxs]
    left_right_pairs = get_left_right_pairs(hapmap_cm, annotation_cm)

    ldscores = []
    for i, (left_idx, right_idx) in enumerate(left_right_pairs):
        variances = annotation_matrix[left_idx:right_idx + 1, :].dot(hapmap_matrix[i, :]) / people
        ldscores.append(((people - 1) * variances.dot(variances) - len(variances)) / (people - 2))

    rs_positions = get_rs_positions(ancestry, chromosome, hapmap_set)
    allele_frequency = np.sum(annotation_matrix, 1) / 2 / people
    minor_allele_frequency = np.minimum(allele_frequency, np.ones(annotation_matrix.shape[0]) - allele_frequency)
    m = annotation_matrix.shape[0]
    m_5_50 = np.sum(minor_allele_frequency > 0.05)

    write_ldscore(dataset, ancestry, chromosome, ldscores, rs_positions)
    write_m(dataset, ancestry, chromosome, m, m_5_50)


def zip_up_data(dataset: str, ancestry: str) -> None:
    os.makedirs(f'data/zip/', exist_ok=True)
    annotations = f'data/annotation/{dataset}.{ancestry}.*.annot.gz'
    ldscores = f'data/ldscore/{dataset}.{ancestry}.*.l2.*'
    subprocess.check_call(f'zip -j data/zip/ldscore.{dataset}.{ancestry}.zip {annotations} {ldscores}', shell=True)


def upload(username: str, dataset: str) -> None:
    path = f'{s3_bucket}/userdata/{username}/annotation/{dataset}/ldsc/ldscore/'
    subprocess.check_call(f'aws s3 cp data/zip/ {path} --recursive', shell=True)


def clean_up() -> None:
    if os.path.exists('data'):
        shutil.rmtree('data')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--username', default=None, required=True, type=str)
    parser.add_argument('--dataset', default=None, required=True, type=str)
    args = parser.parse_args()

    for ancestry in ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']:
        download(args.username, args.dataset, ancestry)
        for chromosome in range(1, 23):
            run(args.dataset, ancestry, chromosome)
        zip_up_data(args.dataset, ancestry)
    upload(args.username, args.dataset)
    clean_up()


if __name__ == '__main__':
    main()
