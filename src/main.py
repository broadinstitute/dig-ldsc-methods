import argparse
import os
import subprocess
import tempfile

s3_bucket = os.environ['S3_BUCKET']
input_path = os.environ['INPUT_PATH']


def download_sumstats(username: str, dataset: str, tmp_dir: str) -> None:
    path = f'{s3_bucket}/userdata/{username}/genetic/{dataset}/raw/'
    subprocess.check_call(f'aws s3 cp {path} {tmp_dir}/dataset/ --recursive', shell=True)


def download_regression(username: str, dataset: str, tmp_dir: str) -> None:
    file = f'{s3_bucket}/userdata/{username}/genetic/{dataset}/ldsc/sumstats/{dataset}.sumstats.gz'
    metadata = f'{s3_bucket}/userdata/{username}/genetic/{dataset}/ldsc/sumstats/metadata'
    subprocess.check_call(f'aws s3 cp {file} {tmp_dir}/dataset/input.sumstats.gz', shell=True)
    subprocess.check_call(f'aws s3 cp {metadata} {tmp_dir}/dataset/', shell=True)


def upload_sumstats(username: str, dataset: str, tmp_dir: str) -> None:
    path = f'{s3_bucket}/userdata/{username}/genetic/{dataset}/ldsc/sumstats/'
    subprocess.check_call(f'aws s3 cp {tmp_dir}/output/ {path} --recursive', shell=True)


def upload_regression(username: str, dataset: str, tmp_dir: str) -> None:
    path = f'{s3_bucket}/userdata/{username}/genetic/{dataset}/ldsc/s-ldsc/'
    subprocess.check_call(f'aws s3 cp {tmp_dir}/output/ {path} --recursive', shell=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--username', default=None, required=True, type=str)
    parser.add_argument('--dataset', default=None, required=True, type=str)
    parser.add_argument('--method', default=None, required=True, type=str)
    args = parser.parse_args()

    tmp = tempfile.TemporaryDirectory()
    if args.method == 'sumstats':
        download_sumstats(args.username, args.dataset, tmp.name)
        subprocess.check_call(['python3.8', 'main.py', f'--dir={tmp.name}'], cwd='ldsc/sumstats/')
        upload_sumstats(args.username, args.dataset, tmp.name)
    if args.method == 'regression':
        download_regression(args.username, args.dataset, tmp.name)
        subprocess.check_call(['python3.8', 'main.py', f'--dir={tmp.name}'], cwd='ldsc/regression/')
        upload_regression(args.username, args.dataset, tmp.name)
    tmp.cleanup()


if __name__ == '__main__':
    main()
