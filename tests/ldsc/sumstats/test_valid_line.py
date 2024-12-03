from typing import Dict
from src.ldsc.sumstats.main import valid_line, var_id_columns


def get_valid_line_data() -> (Dict, Dict):
    line = {'F1': '1', 'F2': 0, 'F3': 'A', 'F4': 'T', 'F5': 1E-5, 'F6': 2, 'F7': 0.8, 'F8': 1000.0}
    col_map = {'chromosome': 'F1', 'position': 'F2', 'reference': 'F3', 'alt': 'F4',
               'pValue': 'F5', 'beta': 'F6', 'oddsRatio': 'F7', 'n': 'F8'}
    return line, col_map


def test_valid_line() -> None:
    line, col_map = get_valid_line_data()
    assert valid_line(line, col_map, None)


def test_invalid_line_missing_var_id_columns() -> None:
    for column in var_id_columns:
        line, col_map = get_valid_line_data()
        col_map[column] = 'fake'
        assert valid_line(line, col_map, None) is False


def test_invalid_line_missing_pValue() -> None:
    line, col_map = get_valid_line_data()
    col_map['pValue'] = 'fake'
    assert valid_line(line, col_map, None) is False

    line, col_map = get_valid_line_data()
    line[col_map['pValue']] = 0
    assert valid_line(line, col_map, None) is False

    line, col_map = get_valid_line_data()
    line[col_map['pValue']] = -1E-5
    assert valid_line(line, col_map, None) is False

    line, col_map = get_valid_line_data()
    line[col_map['pValue']] = 1 + 1E-15
    assert valid_line(line, col_map, None) is False


def test_invalid_line_missing_beta_and_odds_ratio() -> None:
    line, col_map = get_valid_line_data()
    col_map['beta'] = 'fake'
    col_map['oddsRatio'] = 'fake'
    assert valid_line(line, col_map, None) is False

    line, col_map = get_valid_line_data()
    col_map.pop('beta')
    col_map['oddsRatio'] = 'fake'
    assert valid_line(line, col_map, None) is False

    line, col_map = get_valid_line_data()
    col_map['beta'] = 'fake'
    col_map.pop('oddsRatio')
    assert valid_line(line, col_map, None) is False

    line, col_map = get_valid_line_data()
    col_map.pop('beta')
    col_map.pop('oddsRatio')
    assert valid_line(line, col_map, None) is False

    line, col_map = get_valid_line_data()
    col_map['beta'] = 'fake'
    assert valid_line(line, col_map, None)

    line, col_map = get_valid_line_data()
    col_map.pop('beta')
    assert valid_line(line, col_map, None)

    line, col_map = get_valid_line_data()
    col_map['oddsRatio'] = 'fake'
    assert valid_line(line, col_map, None)

    line, col_map = get_valid_line_data()
    col_map.pop('oddsRatio')
    assert valid_line(line, col_map, None)


def test_invalid_line_missing_n() -> None:
    line, col_map = get_valid_line_data()
    col_map['n'] = 'fake'
    assert valid_line(line, col_map, None) is False

    line, col_map = get_valid_line_data()
    col_map['n'] = 'fake'
    assert valid_line(line, col_map, 1000.0)

    line, col_map = get_valid_line_data()
    col_map.pop('n')
    assert valid_line(line, col_map, None) is False

    line, col_map = get_valid_line_data()
    col_map.pop('n')
    assert valid_line(line, col_map, 1000.0)
