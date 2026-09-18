import gzip
from tempfile import TemporaryDirectory

from src.ldsc.sldsc.make_sumstats import stream_to_data

HEADER = ['chromosome', 'position', 'reference', 'alt', 'pValue', 'beta', 'n']
ROWS = [['1', '1', 'A', 'B', '1E-5', '2.0', '10000']]
METADATA = {'col_map': {h: h for h in HEADER}, 'separator': '\t'}
VAR_TO_RS = {'1:1:A:B': 'rs_1'}


def write(path: str, opener) -> None:
    with opener(path, 'wt') as f:
        f.write('\t'.join(HEADER) + '\n')
        for row in ROWS:
            f.write('\t'.join(row) + '\n')


def test_gzipped_file_is_read() -> None:
    tmp = TemporaryDirectory()
    write(f'{tmp.name}/test.tsv.gz', gzip.open)
    out, count = stream_to_data(f'{tmp.name}/test.tsv.gz', VAR_TO_RS, {}, METADATA)
    assert count['translated'] == 1
    assert out[0][0] == 'rs_1'
    tmp.cleanup()


def test_uncompressed_file_is_read() -> None:
    tmp = TemporaryDirectory()
    write(f'{tmp.name}/test.tsv', open)
    out, count = stream_to_data(f'{tmp.name}/test.tsv', VAR_TO_RS, {}, METADATA)
    assert count['translated'] == 1
    assert out[0][0] == 'rs_1'
    tmp.cleanup()
