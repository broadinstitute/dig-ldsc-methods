import glob
import gzip
from multiprocessing import Pool
import numpy as np
from numpy import typing as npt
import os
import re


def frq_path(ancestry: str, chromosome: int) -> str:
    return f'frq/{ancestry}/chr.{chromosome}.frq'


def baseline_annot_path(ancestry: str, chromosome: int) -> str:
    return f'baseline/{ancestry}/baselineLD.{chromosome}.annot.gz'


def tissue_annot_path(tissue: str, ancestry: str, chromosome: int) -> str:
    return f'tissue/{ancestry}/{tissue}/tissueLD.{chromosome}.annot.gz'


def get_baseline_shape(ancestry: str) -> int:
    baseline_file = baseline_annot_path(ancestry, 1)
    with gzip.open(baseline_file, 'rt') as f:
        return len(f.readline().strip().split('\t')) - 4


def get_tissue_shape(tissue: str, ancestry: str) -> int:
    tissue_file = tissue_annot_path(tissue, ancestry, 1)
    with gzip.open(tissue_file, 'rt') as f:
        return len(f.readline().strip().split('\t'))


def run_baseline_overlap(ancestry: str) -> None:
    shape = get_baseline_shape(ancestry)
    output = np.zeros((shape, shape))
    for chromosome in range(1, 23):
        print(chromosome)
        with gzip.open(baseline_annot_path(ancestry, chromosome), 'rt') as file, \
                open(frq_path(ancestry, chromosome), 'r') as frq_file:
            file.readline(), frq_file.readline()
            for line, frq_line in zip(file, frq_file):
                frq = float(frq_line.strip().split('\t')[4])
                if 0.05 < frq < 0.95:
                    values = line.strip().split('\t')[4:]
                    np_array = np.array([list(map(float, values))])
                    output += np_array.T.dot(np_array)
    save_overlap(f'overlap.baseline.{ancestry}.npy', output)


def run_tissue_overlap(file: str) -> None:
    ancestry, tissue = re.findall('tissue/(.*)/(.*)', file)[0]
    baseline_shape = get_baseline_shape(ancestry)
    tissue_shape = get_tissue_shape(tissue, ancestry)
    baseline_tissue_output = np.zeros((baseline_shape, tissue_shape))
    tissue_output = np.zeros((tissue_shape, tissue_shape))
    for chromosome in range(1, 23):
        print(tissue, ancestry, chromosome)
        with gzip.open(baseline_annot_path(ancestry, chromosome), 'rt') as baseline_file, \
                gzip.open(tissue_annot_path(tissue, ancestry, chromosome), 'rt') as tissue_file, \
                open(frq_path(ancestry, chromosome), 'r') as frq_file:
            baseline_file.readline(), tissue_file.readline(), frq_file.readline()
            for baseline_line, tissue_line, frq_line in zip(baseline_file, tissue_file, frq_file):
                frq = float(frq_line.strip().split('\t')[4])
                if 0.05 < frq < 0.95:
                    baseline_values = baseline_line.strip().split('\t')[4:]
                    baseline_np_array = np.array([list(map(float, baseline_values))])
                    tissue_values = tissue_line.strip().split('\t')
                    tissue_np_array = np.array([list(map(float, tissue_values))])
                    baseline_tissue_output += baseline_np_array.T.dot(tissue_np_array)
                    tissue_output += tissue_np_array.T.dot(tissue_np_array)
    save_overlap(f'overlap.baseline.{tissue}.{ancestry}.npy', baseline_tissue_output)
    save_overlap(f'overlap.{tissue}.{ancestry}.npy', tissue_output)



def save_overlap(file: str, overlap: npt.NDArray) -> None:
    if not os.path.exists('datasets'):
        os.mkdir('datasets')
    np.save(f'datasets/{file}', overlap)


def main():
    for ancestry in ['EAS', 'SAS']:
        run_baseline_overlap(ancestry)
        with Pool(10) as p:
            p.map(run_tissue_overlap, glob.glob(f'tissue/{ancestry}/*'))


if __name__ == '__main__':
    main()
