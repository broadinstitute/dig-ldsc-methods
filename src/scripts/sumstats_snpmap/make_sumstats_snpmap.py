import os
import re
import subprocess


group_to_num = {('A', 'C'): 0, ('C', 'A'): 0, ('T', 'G'): 0, ('G', 'T'): 0,
                ('A', 'G'): 1, ('G', 'A'): 1, ('T', 'C'): 1, ('C', 'T'): 1}
num_to_group = [[('A', 'C'), ('C', 'A'), ('T', 'G'), ('G', 'T')], [('A', 'G'), ('G', 'A'), ('T', 'C'), ('C', 'T')]]


def get_hapmap_set():
    out = set()
    with open('data/hapmap/w_hm3.snplist', 'r') as f:
        f.readline()
        for line in f:
            split_line = line.strip().split('\t')
            out |= {split_line[0]}
    return out


def get_g1000_map(ancestry):
    out = {}
    rs_ids = []
    for chromosome in range(1, 23):
        with open(f'data/g1000/{ancestry}/chr{chromosome}.bim', 'r') as f:
            for line in f:
                split_line = line.strip().split('\t')
                rs_id = split_line[1]
                ref_alt = (split_line[5], split_line[4])
                if ref_alt in group_to_num:
                    ref_alts = num_to_group[group_to_num[ref_alt]]
                    out[rs_id] = [f'{split_line[0]}:{split_line[3]}:{ref}:{alt}' for ref, alt in ref_alts]
                    rs_ids.append(rs_id)
    return out, rs_ids


def get_liftover_g1000_map(g1000_map):
    with open('g1000_map.bed', 'w') as f:
        for rs_id, var_ids in g1000_map.items():
            for var_id in var_ids:
                chromosome, position = re.search(r'([^:]*):([^:]*):[^:]*:[^:]*', var_id).groups()
                f.write(f'{chromosome}\t{position}\t{position}\t{rs_id}:{var_id}\n')
    subprocess.check_call('./liftover/liftOver.macOS.arm64 g1000_map.bed liftover/b37toHg38.over.chain hg38.bed out.unmapped', shell=True)

    hg38_g1000_map = {}
    with open('hg38.bed', 'r') as f:
        for line in f:
            chr_hg38_chromosome, hg38_position, _, original_data = line.strip().split('\t')
            hg38_chromosome = chr_hg38_chromosome.replace('chr', '')
            rs_id, chromosome, position, ref, alt = re.search(r'([^:]*):([^:]*):([^:]*):([^:]*):([^:]*)', original_data).groups()
            if hg38_chromosome == chromosome:
                if rs_id not in hg38_g1000_map:
                    hg38_g1000_map[rs_id] = []
                hg38_g1000_map[rs_id].append(f'{chromosome}:{hg38_position}:{ref}:{alt}')
    os.remove('hg38.bed')
    os.remove('g1000_map.bed')
    os.remove('out.unmapped')
    return hg38_g1000_map


def make_snpmap(genome_build, ancestry, hapmap_set, g1000_map, rs_ids):
    os.makedirs(f'data/snpmap/', exist_ok=True)
    with open(f'data/snpmap/sumstats.{genome_build}.{ancestry}.snpmap', 'w') as f:
        for rs_id in rs_ids:
            if rs_id in hapmap_set and rs_id in g1000_map:
                for var_id in g1000_map[rs_id]:
                    f.write(f'{var_id}\t{rs_id}\n')


def main():
    hapmap_set = get_hapmap_set()
    for ancestry in ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']:
        g1000_map, rs_ids = get_g1000_map(ancestry)
        hg38_g1000_map = get_liftover_g1000_map(g1000_map)
        make_snpmap('GRCh37', ancestry, hapmap_set, g1000_map, rs_ids)
        make_snpmap('GRCh38', ancestry, hapmap_set, hg38_g1000_map, rs_ids)


if __name__ == '__main__':
    main()
