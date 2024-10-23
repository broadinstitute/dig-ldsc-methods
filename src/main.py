import argparse
import os
import shutil
import subprocess

s3_bucket = os.environ['S3_BUCKET']


def download_sumstats(username: str, dataset: str) -> None:
    path = f'{s3_bucket}/userdata/{username}/genetic/{dataset}/raw/'
    subprocess.check_call(f'aws s3 cp {path} ldsc/sumstats/dataset/ --recursive', shell=True)


def download_regression(username: str, dataset: str) -> None:
    file = f'{s3_bucket}/userdata/{username}/genetic/{dataset}/ldsc/sumstats/{dataset}.sumstats.gz'
    metadata = f'{s3_bucket}/userdata/{username}/genetic/{dataset}/raw/metadata'
    subprocess.check_call(f'aws s3 cp {file} ldsc/regression/dataset/input.sumstats.gz', shell=True)
    subprocess.check_call(f'aws s3 cp {metadata} ldsc/regression/dataset/', shell=True)


def upload_sumstats(username: str, dataset: str):
    file = 'ldsc/sumstats/sumstats/output.sumstats.gz'
    path = f'{s3_bucket}/userdata/{username}/genetic/{dataset}/ldsc/sumstats/{dataset}.sumstats.gz'
    subprocess.check_call(f'aws s3 cp {file} {path}', shell=True)


def upload_regression(username: str, dataset: str) -> None:
    path = f'{s3_bucket}/userdata/{username}/genetic/{dataset}/ldsc/s-ldsc/'
    files = f'ldsc/regression/regression/'
    subprocess.check_call(f'aws s3 cp {files} {path} --recursive', shell=True)


def clean_up_sumstats():
    for directory in ['ldsc/sumstats/dataset', 'ldsc/sumstats/sumstats']:
        if os.path.exists(directory):
            shutil.rmtree(directory)


def clean_up_regression():
    for directory in ['ldsc/regression/dataset', 'ldsc/regression/regression']:
        if os.path.exists(directory):
            shutil.rmtree(directory)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--username', default=None, required=True, type=str)
    parser.add_argument('--dataset', default=None, required=True, type=str)
    parser.add_argument('--method', default=None, required=True, type=str)
    args = parser.parse_args()

    if args.method == 'sumstats':
        download_sumstats(args.username, args.dataset)
        subprocess.check_call(['python3.8', 'main.py'], cwd='ldsc/sumstats/')
        upload_sumstats(args.username, args.dataset)
        clean_up_sumstats()
    if args.method == 'regression':
        download_regression(args.username, args.dataset)
        subprocess.check_call(['python3.8', 'main.py'], cwd='ldsc/regression/')
        upload_regression(args.username, args.dataset)
        clean_up_regression()


if __name__ == '__main__':
    main()
