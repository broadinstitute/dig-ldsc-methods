import argparse
import json
from typing import Dict

import sldsc, annot_sldsc, make_annot, make_ld, make_sumstats


def get_metadata(data_path: str) -> Dict:
    with open(f'{data_path}/raw/metadata', 'r') as f:
        metadata = json.load(f)
    return metadata


def save_metadata(data_path: str, method: str, metadata: Dict) -> None:
    with open(f'{data_path}/sldsc/{method}/metadata', 'w') as f:
        json.dump(metadata, f)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir', default=None, required=True, type=str)
    parser.add_argument('--method', default=None, required=True, type=str)
    args = parser.parse_args()

    metadata = get_metadata(args.dir)

    if args.method == 'sldsc':
        metadata = make_sumstats.sumstats(args.dir, metadata)
        metadata = sldsc.sldsc(args.dir, metadata)
    elif args.method == 'annot-sldsc':
        make_annot.annotation(args.dir, metadata)
        make_ld.ld(args.dir, metadata)
        metadata = annot_sldsc.annot_sldsc(args.dir, metadata)
    save_metadata(args.dir, args.method, metadata)


if __name__ == '__main__':
    main()
