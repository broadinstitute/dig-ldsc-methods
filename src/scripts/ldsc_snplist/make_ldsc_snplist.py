import os


group_to_num = {('A', 'C'): 0, ('C', 'A'): 0, ('T', 'G'): 0, ('G', 'T'): 0,
                ('A', 'G'): 1, ('G', 'A'): 1, ('T', 'C'): 1, ('C', 'T'): 1}
num_to_group = [[('A', 'C'), ('C', 'A'), ('T', 'G'), ('G', 'T')], [('A', 'G'), ('G', 'A'), ('T', 'C'), ('C', 'T')]]


def get_hapmap_set():
    out = set()
    with open('data/snps/w_hm3.snplist', 'r') as f:
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


def make_snplist(ancestry, hapmap_set, g1000_map, rs_ids):
    if not os.path.exists('data/snplist'):
        os.mkdir('data/snplist')
    with open(f'data/snplist/ldsc.{ancestry}.snplist', 'w') as f:
        for rs_id in rs_ids:
            if rs_id in hapmap_set and rs_id in g1000_map:
                for var_id in g1000_map[rs_id]:
                    f.write(f'{var_id}\t{rs_id}\n')


def main():
    hapmap_set = get_hapmap_set()
    for ancestry in ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']:
        g1000_map, rs_ids = get_g1000_map(ancestry)
        make_snplist(ancestry, hapmap_set, g1000_map, rs_ids)


if __name__ == '__main__':
    main()
