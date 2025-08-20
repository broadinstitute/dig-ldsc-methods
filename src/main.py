import argparse
import os
import subprocess
import tempfile

s3_bucket = os.environ['S3_BUCKET']
input_path = os.environ['INPUT_PATH']


def get_input_output_method(method: str) -> (str, str, str):
    return {
        'sumstats': ('genetic', 'sldsc', 'sumstats'),
        'sldsc': ('genetic', 'sldsc', 'sldsc'),
        'magma-sumstats': ('genetic', 'magma', 'sumstats'),
        'magma-genes': ('genetic', 'magma', 'genes')
    }[method]


def download_raw(username: str, dataset: str, tmp_dir: str) -> None:
    path = f'{s3_bucket}/userdata/{username}/genetic/{dataset}/raw/'
    subprocess.check_call(f'aws s3 cp "{path}" {tmp_dir}/raw/ --recursive', shell=True)


def download_sumstats(username: str, dataset: str, tmp_dir: str) -> None:
    path = f'{s3_bucket}/userdata/{username}/genetic/{dataset}/sldsc/sumstats/'
    subprocess.check_call(f'aws s3 cp "{path}" {tmp_dir}/sumstats/ --recursive', shell=True)


def download_magma_sumstats(username: str, dataset: str, tmp_dir: str) -> None:
    path = f'{s3_bucket}/userdata/{username}/genetic/{dataset}/magma/sumstats/'
    subprocess.check_call(f'aws s3 cp "{path}" {tmp_dir}/sumstats/ --recursive', shell=True)


def upload(username: str, dataset: str, tmp_dir: str, method: str) -> None:
    input_type, output_group, output_method = get_input_output_method(method)
    path = f'{s3_bucket}/userdata/{username}/{input_type}/{dataset}/{output_group}/{output_method}'
    subprocess.check_call(f'aws s3 cp {tmp_dir}/{output_method}/ "{path}" --recursive', shell=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--username', default=None, required=True, type=str)
    parser.add_argument('--dataset', default=None, required=True, type=str)
    parser.add_argument('--method', default=None, required=True, type=str)
    args = parser.parse_args()

    tmp = tempfile.TemporaryDirectory()
    print(f'Using temporary directory {tmp.name}')
    if args.method == 'sumstats':
        download_raw(args.username, args.dataset, tmp.name)
        subprocess.check_call(['python3.8', 'main.py', f'--dir={tmp.name}'], cwd='ldsc/sumstats/')
        upload(args.username, args.dataset, tmp.name, args.method)
    if args.method == 'sldsc':
        download_sumstats(args.username, args.dataset, tmp.name)
        subprocess.check_call(['python3.8', 'main.py', f'--dir={tmp.name}'], cwd='ldsc/sldsc/')
        upload(args.username, args.dataset, tmp.name, args.method)
    if args.method == 'magma-sumstats':
        download_raw(args.username, args.dataset, tmp.name)
        subprocess.check_call(['python3.8', 'main.py', f'--dir={tmp.name}'], cwd='magma/sumstats/')
        upload(args.username, args.dataset, tmp.name, args.method)
    if args.method == 'magma-genes':
        download_magma_sumstats(args.username, args.dataset, tmp.name)
        subprocess.check_call(['python3.8', 'main.py', f'--dir={tmp.name}'], cwd='magma/genes/')
        upload(args.username, args.dataset, tmp.name, args.method)


if __name__ == '__main__':
    main()
