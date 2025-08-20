import os
import re
import subprocess


def get_rs_map():
    pos = {}
    with open(f'files/dbSNP_common_GRCh37.csv', 'r') as f:
        _ = f.readline()
        for line in f:
            split_line = line.strip().split('\t')
            rs_id = split_line[0]
            chromosome, position, ref, alt = split_line[1].split(':')
            if rs_id not in pos:
                pos[chromosome, int(position)] = rs_id
    return pos


def get_liftover_map(pos):
    with open('pos_map.bed', 'w') as f:
        for (chromosome, position), rs_id in pos.items():
            f.write(f'{chromosome}\t{position}\t{position}\t{rs_id}:{chromosome}:{position}\n')
    subprocess.check_call('./liftover/liftOver.macOS.arm64 pos_map.bed liftover/b37toHg38.over.chain hg38.bed out.unmapped', shell=True)

    hg38_pos = {}
    with open('hg38.bed', 'r') as f:
        for line in f:
            chr_hg38_chromosome, hg38_position, _, original_data = line.strip().split('\t')
            hg38_chromosome = chr_hg38_chromosome.replace('chr', '')
            rs_id, chromosome, position = re.search(r'([^:]*):([^:]*):([^:]*)', original_data).groups()
            if hg38_chromosome == chromosome:
                hg38_pos[(chromosome, int(hg38_position))] = rs_id
    os.remove('hg38.bed')
    os.remove('pos_map.bed')
    os.remove('out.unmapped')
    return hg38_pos


def make_snpmap(genome_build, pos):
    os.makedirs(f'data/snpmap/', exist_ok=True)
    with open(f'data/snpmap/sumstats.{genome_build}.snpmap', 'w') as f:
        for chromosome, position in sorted(pos):
            rs_id = pos[(chromosome, position)]
            f.write(f'{chromosome}\t{position}\t{rs_id}\n')


def main():
    pos = get_rs_map()
    hg38_pos = get_liftover_map(pos)
    make_snpmap('GRCh37', pos)
    make_snpmap('GRCh38', hg38_pos)


if __name__ == '__main__':
    main()
