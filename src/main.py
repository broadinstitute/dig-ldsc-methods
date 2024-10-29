import argparse
import os
import subprocess
import tempfile

s3_bucket = os.environ['S3_BUCKET']
input_path = os.environ['INPUT_PATH']


def download_sumstats(username: str, dataset: str, tmp_dir: str) -> None:
    path = f'{s3_bucket}/userdata/{username}/genetic/{dataset}/raw/'
    subprocess.check_call(f'aws s3 cp {path} {tmp_dir}/raw/ --recursive', shell=True)


def download_sldsc(username: str, dataset: str, tmp_dir: str) -> None:
    path = f'{s3_bucket}/userdata/{username}/genetic/{dataset}/sldsc/sumstats/'
    subprocess.check_call(f'aws s3 cp {path} {tmp_dir}/sumstats/ --recursive', shell=True)


def upload_sumstats(username: str, dataset: str, tmp_dir: str) -> None:
    path = f'{s3_bucket}/userdata/{username}/genetic/{dataset}/sldsc/sumstats'
    subprocess.check_call(f'aws s3 cp {tmp_dir}/sumstats/ {path} --recursive', shell=True)


def upload_sldsc(username: str, dataset: str, tmp_dir: str) -> None:
    path = f'{s3_bucket}/userdata/{username}/genetic/{dataset}/sldsc/sldsc/'
    subprocess.check_call(f'aws s3 cp {tmp_dir}/sldsc/ {path} --recursive', shell=True)


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
    if args.method == 'sldsc':
        download_sldsc(args.username, args.dataset, tmp.name)
        subprocess.check_call(['python3.8', 'main.py', f'--dir={tmp.name}'], cwd='ldsc/sldsc/')
        upload_sldsc(args.username, args.dataset, tmp.name)
    tmp.cleanup()


if __name__ == '__main__':
    main()
