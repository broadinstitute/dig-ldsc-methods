import argparse
import json
from typing import Dict

import annotation, ld


def get_metadata(data_path: str) -> Dict:
    with open(f'{data_path}/raw/metadata', 'r') as f:
        metadata = json.load(f)
    return metadata


def save_metadata(data_path: str, metadata: Dict) -> None:
    with open(f'{data_path}/sldsc/ld/metadata', 'w') as f:
        json.dump(metadata, f)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir', default=None, required=True, type=str)
    data_path = parser.parse_args().dir

    metadata = get_metadata(data_path)
    ancestry = metadata['ancestry']

    annotation.annotation(data_path, ancestry)
    ld.ld(data_path, ancestry)
    save_metadata(data_path, metadata)


if __name__ == '__main__':
    main()
