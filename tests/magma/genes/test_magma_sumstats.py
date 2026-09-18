import gzip
from tempfile import TemporaryDirectory

from src.magma.genes.sumstats import stream_to_data

HEADER = ['chromosome', 'position', 'pValue', 'n']
ROWS = [['1', '1', '1E-5', '10000']]
METADATA = {'col_map': {h: h for h in HEADER}, 'separator': '\t'}
RS_MAP = {('1', '1'): 'rs_1'}


def write(path: str, opener) -> None:
    with opener(path, 'wt') as f:
        f.write('\t'.join(HEADER) + '\n')
        for row in ROWS:
            f.write('\t'.join(row) + '\n')


def test_gzipped_file_is_read() -> None:
    tmp = TemporaryDirectory()
    write(f'{tmp.name}/test.tsv.gz', gzip.open)
    out, count = stream_to_data(f'{tmp.name}/test.tsv.gz', RS_MAP, METADATA)
    assert count['final'] == 1
    assert out[0] == ('rs_1', 1E-5, 10000.0)
    tmp.cleanup()


def test_uncompressed_file_is_read() -> None:
    tmp = TemporaryDirectory()
    write(f'{tmp.name}/test.tsv', open)
    out, count = stream_to_data(f'{tmp.name}/test.tsv', RS_MAP, METADATA)
    assert count['final'] == 1
    assert out[0] == ('rs_1', 1E-5, 10000.0)
    tmp.cleanup()
