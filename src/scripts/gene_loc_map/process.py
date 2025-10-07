

def get_gene_loc_map():
    gene_loc_map = {}
    with open('files/NCBI37.3.plink.gene.loc', 'r') as f:
        for line in f:
            split_line = line.strip().split('\t')
            chromosome, start, end = split_line[1:4]
            gene = split_line[-1]
            gene_loc_map[gene] = (chromosome, int(start), int(end))
    return gene_loc_map


def add_synonyms(gene_loc_map):
    with open('files/gencode.gene.map', 'r') as f:
        for line in f:
            ensembl, gene = line.strip().split('\t')
            if ensembl not in gene_loc_map and gene in gene_loc_map:
                gene_loc_map[ensembl] = gene_loc_map[gene]
            if ensembl in gene_loc_map and gene not in gene_loc_map:
                gene_loc_map[gene] = gene_loc_map[ensembl]
    return gene_loc_map


def write_file(gene_loc_map):
    with open('processed/gene.loc', 'w') as f:
        for gene, (chromosome, start, end) in gene_loc_map.items():
            f.write('{}\t{}\t{}\t{}\n'.format(gene, chromosome, start, end))


def main():
    gene_loc_map = get_gene_loc_map()
    gene_loc_map = add_synonyms(gene_loc_map)
    write_file(gene_loc_map)


if __name__ == '__main__':
    main()
