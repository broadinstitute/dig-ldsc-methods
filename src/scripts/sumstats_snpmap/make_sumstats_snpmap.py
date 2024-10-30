import os
import re
import subprocess

group_to_num = {('A', 'C'): 0, ('T', 'G'): 0,
                ('C', 'A'): 1, ('G', 'T'): 1,
                ('A', 'G'): 2, ('T', 'C'): 2,
                ('G', 'A'): 3, ('C', 'T'): 3}
num_to_group = [[('A', 'C'), ('T', 'G')], [('C', 'A'), ('G', 'T')], [('A', 'G'), ('T', 'C')], [('G', 'A'),  ('C', 'T')]]
num_to_flip = [[('C', 'A'),  ('G', 'T')], [('A', 'C'), ('T', 'G')], [('G', 'A'),  ('C', 'T')], [('A', 'G'), ('T', 'C')]]


def get_hapmap_set():
    out = set()
    with open('data/hapmap/w_hm3.snplist', 'r') as f:
        f.readline()
        for line in f:
            split_line = line.strip().split('\t')
            out |= {split_line[0]}
    return out


def get_g1000_maps(ancestry):
    pos = {}
    ref_alts = {}
    flipped = {}
    rs_ids = []
    for chromosome in range(1, 23):
        with open(f'data/g1000/{ancestry}/chr{chromosome}.bim', 'r') as f:
            for line in f:
                split_line = line.strip().split('\t')
                rs_id = split_line[1]
                ref_alt = (split_line[5], split_line[4])
                if ref_alt in group_to_num:
                    num = group_to_num[ref_alt]
                    pos[rs_id] = f'{split_line[0]}:{split_line[3]}'
                    ref_alts[rs_id] = [f'{ref}:{alt}' for ref, alt in num_to_group[num]]
                    flipped[rs_id] = [f'{ref}:{alt}' for ref, alt in num_to_flip[num]]
                    rs_ids.append(rs_id)
    return pos, ref_alts, flipped, rs_ids


def get_liftover_g1000_map(g1000_pos):
    with open('g1000_map.bed', 'w') as f:
        for rs_id, pos in g1000_pos.items():
            chromosome, position = pos.split(':')
            f.write(f'{chromosome}\t{position}\t{position}\t{rs_id}:{pos}\n')
    subprocess.check_call('./liftover/liftOver.macOS.arm64 g1000_map.bed liftover/b37toHg38.over.chain hg38.bed out.unmapped', shell=True)

    hg38_g1000_pos = {}
    with open('hg38.bed', 'r') as f:
        for line in f:
            chr_hg38_chromosome, hg38_position, _, original_data = line.strip().split('\t')
            hg38_chromosome = chr_hg38_chromosome.replace('chr', '')
            rs_id, chromosome, position = re.search(r'([^:]*):([^:]*):([^:]*)', original_data).groups()
            if hg38_chromosome == chromosome:
                hg38_g1000_pos[rs_id] = f'{chromosome}:{hg38_position}'
    os.remove('hg38.bed')
    os.remove('g1000_map.bed')
    os.remove('out.unmapped')
    return hg38_g1000_pos


def make_snpmap(build_type, genome_build, ancestry, hapmap_set, g1000_pos, g1000_map, rs_ids):
    os.makedirs(f'data/snpmap/', exist_ok=True)
    with open(f'data/snpmap/sumstats.{build_type}.{genome_build}.{ancestry}.snpmap', 'w') as f:
        for rs_id in rs_ids:
            if rs_id in hapmap_set and rs_id in g1000_pos and rs_id in g1000_map:
                for var_id in g1000_map[rs_id]:
                    f.write(f'{g1000_pos[rs_id]}:{var_id}\t{rs_id}\n')


def main():
    hapmap_set = get_hapmap_set()
    for ancestry in ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']:
        g1000_pos, g1000_map, g1000_flipped, rs_ids = get_g1000_maps(ancestry)
        hg38_g1000_pos = get_liftover_g1000_map(g1000_pos)
        make_snpmap('standard','GRCh37', ancestry, hapmap_set, g1000_pos, g1000_map, rs_ids)
        make_snpmap('standard','GRCh38', ancestry, hapmap_set, hg38_g1000_pos, g1000_map, rs_ids)
        make_snpmap('flipped','GRCh37', ancestry, hapmap_set, g1000_pos, g1000_flipped, rs_ids)
        make_snpmap('flipped','GRCh38', ancestry, hapmap_set, hg38_g1000_pos, g1000_flipped, rs_ids)


if __name__ == '__main__':
    main()
