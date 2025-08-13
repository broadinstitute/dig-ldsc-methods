import os
import re
import subprocess


conversion = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
def compliment(bases):
    return ''.join([conversion[base] for base in bases])


def get_rs_map():
    pos = {}
    ref_alts = {}
    flipped = {}
    rs_ids = []
    with open(f'files/dbSNP_common_GRCh37.csv', 'r') as f:
        _ = f.readline()
        for line in f:
            split_line = line.strip().split('\t')
            rs_id = split_line[0]
            chromosome, position, ref, alt = split_line[1].split(':')
            if all([a in conversion for a in ref + alt]):
                if rs_id not in pos:
                    pos[rs_id] = f'{chromosome}:{position}'
                    rs_ids.append(rs_id)
                    ref_alts[rs_id] = []
                    flipped[rs_id] = []
                ref_alts[rs_id] += [f'{ref}:{alt}', f'{compliment(ref)}:{compliment(alt)}']
                flipped[rs_id] += [f'{alt}:{ref}', f'{compliment(alt)}:{compliment(ref)}']
            print(len(rs_ids))
    return pos, ref_alts, flipped, rs_ids


def get_liftover_map(pos):
    with open('pos_map.bed', 'w') as f:
        for rs_id, pos in pos.items():
            chromosome, position = pos.split(':')
            f.write(f'{chromosome}\t{position}\t{position}\t{rs_id}:{pos}\n')
    subprocess.check_call('./liftover/liftOver.macOS.arm64 pos_map.bed liftover/b37toHg38.over.chain hg38.bed out.unmapped', shell=True)

    hg38_pos = {}
    with open('hg38.bed', 'r') as f:
        for line in f:
            chr_hg38_chromosome, hg38_position, _, original_data = line.strip().split('\t')
            hg38_chromosome = chr_hg38_chromosome.replace('chr', '')
            rs_id, chromosome, position = re.search(r'([^:]*):([^:]*):([^:]*)', original_data).groups()
            if hg38_chromosome == chromosome:
                hg38_pos[rs_id] = f'{chromosome}:{hg38_position}'
    os.remove('hg38.bed')
    os.remove('pos_map.bed')
    os.remove('out.unmapped')
    return hg38_pos


def make_snpmap(build_type, genome_build, pos, ref_alts, rs_ids):
    os.makedirs(f'data/snpmap/', exist_ok=True)
    with open(f'data/snpmap/sumstats.{build_type}.{genome_build}.snpmap', 'w') as f:
        for rs_id in rs_ids:
            if rs_id in pos and rs_id in ref_alts:
                for var_id in ref_alts[rs_id]:
                    f.write(f'{pos[rs_id]}:{var_id}\t{rs_id}\n')


def main():
    pos, ref_alts, flipped, rs_ids = get_rs_map()
    hg38_pos = get_liftover_map(pos)
    make_snpmap('standard','GRCh37', pos, ref_alts, rs_ids)
    make_snpmap('standard','GRCh38', hg38_pos, ref_alts, rs_ids)
    make_snpmap('flipped','GRCh37', pos, flipped, rs_ids)
    make_snpmap('flipped','GRCh38', hg38_pos, flipped, rs_ids)


if __name__ == '__main__':
    main()
