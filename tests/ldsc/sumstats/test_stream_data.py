import gzip
from tempfile import TemporaryDirectory
from typing import Dict, List
from src.ldsc.sumstats.main import stream_to_data, p_to_z


def get_header():
    return ['chromosome', 'position', 'reference', 'alt', 'pValue', 'beta', 'n']


def get_metadata(sep: str) -> Dict:
    metadata = {}
    metadata['col_map'] = {header: header for header in get_header()}
    metadata['separator'] = sep
    return metadata


def get_var_to_rs_map() -> Dict:
    return {'1:1:A:B': 'rs_1', '1:2:B:A': 'rs_2', '1:3:C:D': 'rs_3'}


def get_flipped_var_to_rs_map() -> Dict:
    return {'1:4:D:C': 'rs_4', '1:5:E:F': 'rs_5', '1:6:F:E': 'rs_6'}


def default_line() -> Dict:
    return {'chromosome': '1', 'position': 1, 'reference': 'A', 'alt': 'B', 'pValue': 1E-5, 'beta': 2.0, 'n': 10000.0}


def get_line(overrides: Dict) -> Dict:
    line = default_line()
    return {header: overrides.get(header, line[header]) for header in get_header()}


def make_file(lines: List, sep: str = ',') -> TemporaryDirectory:
    tmp = TemporaryDirectory()
    headers = get_header()
    with gzip.open(f'{tmp.name}/test.csv.gz', 'wt') as f:
        f.write(f'{sep.join(headers)}\n')
        for line in lines:
            f.write(f'{sep.join([str(line[header]) for header in headers])}\n')
    return tmp


def valid_lines() -> List:
    return [
        get_line({'position': position, 'reference': reference, 'alt': alt, 'pValue':  pValue, 'beta': beta, 'n': n})
        for position, reference, alt, pValue, beta, n in [
            (1, 'A', 'B', 1E-5, 2.0, 10000.0),
            (2, 'B', 'A', 2E-5, 2.1, 11000.0),
            (3, 'C', 'D', 3E-5, 2.2, 12000.0),
            (4, 'D', 'C', 4E-5, 2.3, 13000.0),
            (5, 'E', 'F', 5E-5, 2.4, 14000.0),
            (6, 'F', 'E', 6E-5, 2.5, 15000.0)
        ]
    ]


def test_valid_csv_stream() -> None:
    metadata = get_metadata(',')
    lines = valid_lines()
    tmp = make_file(lines)
    out, count = stream_to_data(f'{tmp.name}/test.csv.gz', get_var_to_rs_map(), get_flipped_var_to_rs_map(), metadata)
    assert count['all'] == 6
    assert count['error'] == 0
    assert count['flipped'] == 3
    assert count['translated'] == 6
    assert out[0] == ('rs_1', p_to_z(1E-5, 2.0), 10000.0)
    assert out[1] == ('rs_2', p_to_z(2E-5, 2.1), 11000.0)
    assert out[2] == ('rs_3', p_to_z(3E-5, 2.2), 12000.0)
    assert out[3] == ('rs_4', p_to_z(4E-5, -2.3), 13000.0)
    assert out[4] == ('rs_5', p_to_z(5E-5, -2.4), 14000.0)
    assert out[5] == ('rs_6', p_to_z(6E-5, -2.5), 15000.0)
    tmp.cleanup()


def test_skipped_csv_stream() -> None:
    metadata = get_metadata(',')
    lines = valid_lines()
    lines[0]['chromosome'] = '2'
    lines[1]['position'] = 20
    lines[4]['reference'] = 'fake'
    lines[5]['alt'] = 'fake'
    tmp = make_file(lines)
    out, count = stream_to_data(f'{tmp.name}/test.csv.gz', get_var_to_rs_map(), get_flipped_var_to_rs_map(), metadata)
    assert count['all'] == 6
    assert count['error'] == 0
    assert count['flipped'] == 1
    assert count['translated'] == 2
    tmp.cleanup()


def test_error_csv_stream() -> None:
    metadata = get_metadata(',')
    lines = valid_lines()
    headers = get_header()
    for header in headers:
        line = get_line({header: ''})
        lines.append(line)
    tmp = make_file(lines)
    out, count = stream_to_data(f'{tmp.name}/test.csv.gz', get_var_to_rs_map(), get_flipped_var_to_rs_map(), metadata)
    assert count['all'] == 6 + len(headers)
    assert count['error'] == 2  # beta and n fail to cast in translation
    assert count['flipped'] == 3
    assert count['translated'] == 6
    tmp.cleanup()


def test_tsv() -> None:
    metadata = get_metadata('\t')
    lines = valid_lines()
    tmp = make_file(lines, '\t')
    out, count = stream_to_data(f'{tmp.name}/test.csv.gz', get_var_to_rs_map(), get_flipped_var_to_rs_map(), metadata)
    assert count['all'] == 6
    assert count['error'] == 0
    assert count['flipped'] == 3
    assert count['translated'] == 6


def test_not_a_csv() -> None:
    metadata = get_metadata(',')
    lines = valid_lines()
    tmp = make_file(lines, '\t')
    out, count = stream_to_data(f'{tmp.name}/test.csv.gz', get_var_to_rs_map(), get_flipped_var_to_rs_map(), metadata)
    assert count['all'] == 6
    assert count['error'] == 0
    assert count['flipped'] == 0
    assert count['translated'] == 0
