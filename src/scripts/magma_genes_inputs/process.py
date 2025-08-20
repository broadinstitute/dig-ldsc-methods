import os
import subprocess


def convert_chromosome(chromosome):
    if chromosome in [str(i) for i in range(1, 23)]:
        return chromosome
    elif chromosome == 'X':
        return '23'
    elif chromosome == 'Y':
        return '24'


def convert_common():
    with open('variants.csv', 'w') as f_out:
        with open('files/dbSNP_common_GRCh37.csv', 'r') as f:
            _ = f.readline()
            for line in f:
                rs_id, var_id = line.strip().split('\t')
                chromosome, position, ref, alt = var_id.split(':')
                converted_chromosome = convert_chromosome(chromosome)
                if converted_chromosome is not None:
                    f_out.write(f'{rs_id}\t{converted_chromosome}\t{position}\n')


# Change to magma-linux if using on linux machine
def assign_genes():
    subprocess.check_call('docker run -v .:/mnt/var/ -it ubuntu:18.04 /mnt/var/files/magma/magma-osx '
                          '--annotate '
                          '--snp-loc /mnt/var/variants.csv '
                          '--gene-loc /mnt/var/files/NCBI37.3.gene.loc '
                          '--out /mnt/var/processed/all', shell=True)


def main():
    convert_common()
    assign_genes()
    os.remove('variants.csv')


if __name__ == '__main__':
    main()
