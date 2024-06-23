from sequence import Sequence
import argparse
import re


def parse_args():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-file', type=str,
                        help='File containing accession numbers')
    parser.add_argument('-type', type=str, help='Type of alignment')
    parser.add_argument('-format', type=str, help='Fileformat')

    if not (parser.parse_args().file) or not (parser.parse_args().type):
        parser.print_help()
        exit()

    return parser.parse_args()


def load_accession_numbers(file, seq: Sequence):
    with open(file) as f:
        for line in f.readlines():
            split = re.split(r'\s+', line)
            seq.add_sequence([split[0], split[1]])


def main():
    parser = parse_args()

    ff = parser.format

    seq = Sequence(parser.type, ff)

    load_accession_numbers(parser.file, seq)

    seq.print_sequences()
    seq.export()
    if parser.type == "nucleotide":
        seq.do_alignment_matrices()
    else:
        seq.do_alignment()


if __name__ == '__main__':
    main()
