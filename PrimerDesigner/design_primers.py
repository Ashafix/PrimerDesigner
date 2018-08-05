import sys
import argparse
from PrimerDesigner import design_primers


def parse_args(args):
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('targetSequence', type=str, help='either a plain sequence or a FASTA file')
    parser.add_argument('primerPairs', type=int, help='number of primer pairs which should be designed')
    return parser.parse_args(args)


def validate_args(args):
    if not os.path.isfile(args.targetSequence):
        print('Could not find input sequence: {}'.format(args.targetSequence), file=sys.stderr)
        return False
    if args.primerPairs < 1:
        print('Primer pairs need to be least 1', file=sys.stderr)
        return False
    return True


if __name__ == '__main__':
    args = parse_args(sys.argv[1:])
    if not validate_args(args):
        sys.exit(1)
    print(args)
    p = design_primers(args.targetSequence, args.primerPairs)
    print(p)

