import argparse
import json
import os
import subprocess
import tempfile
from typing import Dict

s3_bucket = os.environ['S3_BUCKET']
input_path = os.environ['INPUT_PATH']


def get_config(method: str) -> Dict:
    with open('config.json', 'r') as f:
        return json.load(f)['methods'][method]


def download(username: str, dataset: str, tmp_dir: str, config: Dict) -> None:
    path = f'{s3_bucket}/userdata/{username}/{config["dataset_type"]}/{dataset}/{config["download_from"]}/'
    subprocess.check_call(f'aws s3 cp "{path}" {tmp_dir}/{config["download_from"]}/ --recursive', shell=True)


def upload(username: str, dataset: str, tmp_dir: str, config: Dict) -> None:
    path = f'{s3_bucket}/userdata/{username}/{config["dataset_type"]}/{dataset}/{config["upload_to"]}/'
    subprocess.check_call(f'aws s3 cp {tmp_dir}/{config["upload_to"]}/ "{path}" --recursive', shell=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--username', default=None, required=True, type=str)
    parser.add_argument('--dataset', default=None, required=True, type=str)
    parser.add_argument('--method', default=None, required=True, type=str)
    args = parser.parse_args()

    config = get_config(args.method)
    tmp = tempfile.TemporaryDirectory()
    download(args.username, args.dataset, tmp.name, config)
    subprocess.check_call(['python3.8', 'main.py', f'--dir={tmp.name}', f'--method={args.method}'], cwd=config['cwd'])
    upload(args.username, args.dataset, tmp.name, config)


if __name__ == '__main__':
    main()
