import gzip
import os
from tempfile import TemporaryDirectory

from src.pigean.pigean.sumstats import stream_to_sumstats

HEADER = ['chromosome', 'position', 'pValue', 'n']
ROWS = [['1', '1', '1E-5', '10000']]
METADATA = {'col_map': {h: h for h in HEADER}, 'separator': '\t'}


def write(path: str, opener) -> None:
    with opener(path, 'wt') as f:
        f.write('\t'.join(HEADER) + '\n')
        for row in ROWS:
            f.write('\t'.join(row) + '\n')


def read_output(data_path: str) -> list:
    with gzip.open(f'{data_path}/pigean/sumstats/pigean.sumstats.gz', 'rt') as f:
        return f.read().splitlines()


def test_gzipped_file_is_read() -> None:
    tmp = TemporaryDirectory()
    os.makedirs(f'{tmp.name}/raw')
    write(f'{tmp.name}/raw/test.tsv.gz', gzip.open)
    count = stream_to_sumstats(tmp.name, 'test.tsv.gz', METADATA)
    assert count['translated'] == 1
    assert read_output(tmp.name) == ['CHROM\tPOS\tP\tN', '1\t1\t1e-05\t10000.0']
    tmp.cleanup()


def test_uncompressed_file_is_read() -> None:
    tmp = TemporaryDirectory()
    os.makedirs(f'{tmp.name}/raw')
    write(f'{tmp.name}/raw/test.tsv', open)
    count = stream_to_sumstats(tmp.name, 'test.tsv', METADATA)
    assert count['translated'] == 1
    assert read_output(tmp.name) == ['CHROM\tPOS\tP\tN', '1\t1\t1e-05\t10000.0']
    tmp.cleanup()
