import argparse
import json
import os
from typing import Dict

import pigean, sumstats

input_path = os.environ.get('INPUT_PATH')
s3_path = os.environ.get('S3_BUCKET')


def check_envvars():
    assert input_path is not None
    assert s3_path is not None


def get_metadata(data_path: str) -> Dict:
    with open(f'{data_path}/raw/metadata', 'r') as f:
        metadata = json.load(f)
    return metadata


def save_metadata(data_path: str, metadata: Dict) -> None:
    with open(f'{data_path}/pigean/pigean/metadata', 'w') as f:
        json.dump(metadata, f)


def main() -> None:
    check_envvars()
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir', default=None, required=True, type=str)
    parser.add_argument('--method', default=None, required=True, type=str)
    data_path = parser.parse_args().dir

    metadata = get_metadata(data_path)
    metadata['input_type'] = 'sumstats'
    metadata['model'] = 'small'

    metadata = sumstats.main(data_path, metadata)
    if metadata['counts']['translated'] > 0:
        pigean.main(data_path, metadata)
    save_metadata(data_path, metadata)


if __name__ == '__main__':
    main()