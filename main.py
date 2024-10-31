from sequence import Sequence
import argparse

def main(args):
    seq = Sequence(args.database, args.file, args.format, args.file, args.export)
    print(seq)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Align multiples sequences')
    parser.add_argument('-database', type=str, help='Database of alignment', required=True)
    parser.add_argument('-file', type=str, help='File containing accession numbers', required=True)
    parser.add_argument('-format', type=str, help='Fileformat fasta, fa, gb, etc')
    parser.add_argument('-e', '--export', action=argparse.BooleanOptionalAction, help='Export sequences')
    args = parser.parse_args()
    main(args)
