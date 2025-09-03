import time
import gzip
import numpy as np
import os
import struct
import subprocess
from typing import List, Tuple
from numpy import typing as npt

input_path = os.environ.get('INPUT_PATH')
s3_path = os.environ.get('S3_BUCKET')

def check_envvars() -> None:
    assert input_path is not None
    assert s3_path is not None


def g1000_path(ancestry: str, chromosome: str, extension: str) -> str:
    return f'{input_path}/g1000/{ancestry}/chr{chromosome}.{extension}'


def check_g1000(ancestry: str) -> None:
    if not os.path.exists(g1000_path(ancestry, '1', 'bim')):
        cmd = f'./bootstrap/g1000.bootstrap.sh {s3_path} {input_path} {ancestry}'
        subprocess.check_call(cmd, shell=True)


def hapmap_path(chromosome: str) -> str:
    return f'{input_path}/hapmap/hm.{chromosome}.snp'


def check_hapmap() -> None:
    if not os.path.exists(hapmap_path('1')):
        cmd = f'./bootstrap/hapmap.bootstrap.sh {s3_path} {input_path}'
        subprocess.check_call(cmd, shell=True)


def annotation_path(data_path: str, chromosome: str) -> str:
    return f'{data_path}/sldsc/ld/ld.{chromosome}.annot.gz'


def get_bim_data(ancestry: str, chromosome: str) -> Tuple[List[int], List[float], List[Tuple[str, int]]]:
    hm3_out = []
    all_cm = []
    rsids = []
    hm_set = set()
    with open(hapmap_path(chromosome), 'r') as f:
        for line in f.readlines():
            hm_set |= {line.strip()}
    with open(g1000_path(ancestry, chromosome, 'bim'), 'r') as f:
        for idx, line in enumerate(f.readlines()):
            split_line = line.strip().split('\t')
            rsid = split_line[1]
            all_cm.append(float(split_line[2]))
            if rsid in hm_set:
                hm3_out.append(idx)
                rsids.append((rsid, int(split_line[3])))
    return hm3_out, all_cm, rsids


def get_annotation(data_path: str, chromosome: str) -> List[int]:
    with gzip.open(annotation_path(data_path, chromosome), 'rt') as f:
        f.readline() # header
        return [idx for idx, line in enumerate(f.readlines()) if line == '1\n']


def write_output(data_path: str, chromosome: str, l2s: List[float], rsids: List[Tuple[str, int]], M: int, M_5: int) -> None:
    os.makedirs(f'{data_path}/sldsc/ld/', exist_ok=True)
    with gzip.open(f'{data_path}/sldsc/ld/ld.{chromosome}.l2.ldscore.gz', 'w') as f:
        f.write(b'CHR\tSNP\tBP\tL2\n')
        for (rsid, bp), l2 in zip(rsids, l2s):
            f.write(f'{chromosome}\t{rsid}\t{bp}\t{round(l2, 3)}\n'.encode())
    with open(f'{data_path}/sldsc/ld/ld.{chromosome}.l2.M', 'w') as f:
        f.write(f'{M}\n')
    with open(f'{data_path}/sldsc/ld/ld.{chromosome}.l2.M_5_50', 'w') as f:
        f.write(f'{M_5}\n')


def get_dimensions(ancestry: str, chromosome: str) -> Tuple[int, int]:
    with open(g1000_path(ancestry, chromosome, 'fam'), 'r') as f:
        num_people = len(f.readlines())

    extra_columns = (4 - num_people % 4) if num_people % 4 != 0 else 0
    bed_width = num_people + extra_columns

    return bed_width, num_people


decode = {'01': 1, '11': 2, '00': 0, '10':9}
decode_full = {}
for i in range(256):
    z = format(i, '08b')[::-1]
    y = [decode[z[i:i+2]] for i in range(0, len(z), 2)]
    decode_full[i] = y

def get_X(ancestry:str, chromosome: str, bed_width: int, idxs: List[int]) -> npt.NDArray:
    with open(g1000_path(ancestry, chromosome, 'bed'), 'rb') as fh:
        bed = fh.read()
    bf = '<' + 'B' * (bed_width // 4)
    return np.array(
        [[b for a in struct.unpack_from(bf, bed, 3 + idx * (bed_width // 4)) for b in decode_full[a]] for idx in idxs]
    )


def normalize_X(X: npt.NDArray) -> npt.NDArray:
    mean = np.mean(X, 1, keepdims=True)
    std = np.std(X, 1, keepdims=True)
    ones_array = np.ones((1, X.shape[1]))
    return (X - mean.dot(ones_array)) / std.dot(ones_array)


def get_LR(x_cm: List[float], y_cm: List[float], window_cm: float=1.0) -> List[Tuple[int, int]]:
    lr = []
    idx = left_idx = right_idx = 0
    while idx < len(x_cm):
        while right_idx < len(y_cm) and y_cm[right_idx] - x_cm[idx] <= window_cm:
            right_idx += 1
        while x_cm[idx] - y_cm[left_idx] > window_cm:
            left_idx += 1
        lr.append((left_idx, right_idx - 1))
        idx += 1
    return lr


def run_chromosome(data_path: str, ancestry: str, chromosome: str) -> None:
    bed_width, num_people = get_dimensions(ancestry, chromosome)
    hm3_idxs, all_cm, rsids = get_bim_data(ancestry, chromosome)
    annot = get_annotation(data_path, chromosome)

    X = get_X(ancestry, chromosome, bed_width, hm3_idxs)[:, :num_people]
    Y = get_X(ancestry, chromosome, bed_width, annot)[:, :num_people]

    # Filter out fully heterogeneous SNPs
    X_filter = np.sum(X == 1, 1) != num_people
    if sum(X_filter) < X.shape[0]:
        raise Exception('X filter cannot be triggered')

    Y_filter = np.sum(Y == 1, 1) != num_people
    annot = [a for i, a in enumerate(annot) if Y_filter[i]]
    Y = Y[Y_filter, :]

    # This is used for the M_5_50 calculation
    AF = np.sum(Y, 1) / 2 / num_people
    MAF = np.minimum(AF, np.ones(Y.shape[0]) - AF)
    M = Y.shape[0]
    M_5 = np.sum(MAF > 0.05)

    X = normalize_X(X)
    Y = normalize_X(Y)

    # This gets left-right bounds for matrix Y per hapmap3 SNP in the g1000 dataset
    x_cm = np.array(all_cm)[hm3_idxs]
    y_cm = np.array(all_cm)[annot]
    lr = get_LR(x_cm, y_cm)

    l2s = []
    for i, (left_idx, right_idx) in enumerate(lr):
        # LD Score calculation for SNP i
        v = Y[left_idx:right_idx + 1, :].dot(X[i, :]) / num_people
        # L2 adjustment
        l2s.append(((num_people - 1) * v.dot(v) - len(v)) / (num_people - 2))

    write_output(data_path, chromosome, l2s, rsids, M, M_5)


def ld(data_path: str, ancestry: str) -> None:
    check_g1000(ancestry)
    check_hapmap()
    tot_time = 0
    for chromosome in range(1, 22):
        t = time.time()
        run_chromosome(data_path, ancestry, str(chromosome))
        print(f'Chromosome Run Time ({chromosome}): {round(time.time() - t, 3)}')
        tot_time += time.time() - t
    print(f'Full LD Run Time: {tot_time}')
