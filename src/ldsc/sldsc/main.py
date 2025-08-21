import argparse

import sldsc, annot_sldsc

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir', default=None, required=True, type=str)
    parser.add_argument('--method', default=None, required=True, type=str)
    args = parser.parse_args()

    if args.method == 'sldsc':
        sldsc.sldsc(args.dir)
    elif args.method == 'annot-sldsc':
        annot_sldsc.annot_sldsc(args.dir)


if __name__ == '__main__':
    main()