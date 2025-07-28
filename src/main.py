import argparse
import os
import subprocess
import tempfile

s3_bucket = os.environ['S3_BUCKET']
input_path = os.environ['INPUT_PATH']


def get_upload_method(method: str) -> (str, str):
    return {
        'sumstats': ('genetic', 'sldsc/sumstats'),
        'sldsc': ('genetic', 'sldsc/sldsc'),
        'magma-sumstats': ('genetic', 'magma/sumstats'),
        'magma-genes': ('genetic', 'magma/genes')
    }[method]


def get_download_method(method: str) -> (str, str):
    return {
        'sumstats': ('genetic', 'raw'),
        'sldsc': ('genetic', 'sldsc/sumstats'),
        'magma-sumstats': ('genetic', 'raw'),
        'magma-genes': ('genetic', 'magma/sumstats')
    }[method]


def get_cwd(method: str) -> str:
    return {
        'sumstats': 'ldsc/sumstats/',
        'sldsc': 'ldsc/sldsc/',
        'magma-sumstats': 'magma/sumstats/',
        'magma-genes': 'magma/genes/',
        'ld': 'ldsc/ld/',
    }[method]


def download(username: str, dataset: str, tmp_dir: str, method: str) -> None:
    dataset_type, download_type = get_download_method(method)
    path = f'{s3_bucket}/userdata/{username}/{dataset_type}/{dataset}/{download_type}/'
    subprocess.check_call(f'aws s3 cp "{path}" {tmp_dir}/{download_type}/ --recursive', shell=True)


def upload(username: str, dataset: str, tmp_dir: str, method: str) -> None:
    dataset_type, upload_type = get_upload_method(method)
    path = f'{s3_bucket}/userdata/{username}/{dataset_type}/{dataset}/{upload_type}'
    subprocess.check_call(f'aws s3 cp {tmp_dir}/{upload_type}/ "{path}" --recursive', shell=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--username', default=None, required=True, type=str)
    parser.add_argument('--dataset', default=None, required=True, type=str)
    parser.add_argument('--method', default=None, required=True, type=str)
    args = parser.parse_args()

    tmp = tempfile.TemporaryDirectory()
    download(args.username, args.dataset, tmp.name, args.method)
    subprocess.check_call(['python3.8', 'main.py', f'--dir={tmp.name}'], cwd=get_cwd(args.method))
    upload(args.username, args.dataset, tmp.name, args.method)


if __name__ == '__main__':
    main()
