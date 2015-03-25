from phylip import (
    get_seqmapping,
    rename_sequences
)
import argparse

def main():
    args = parse_args()

    # Get the mapping of {Seq#_: origname}
    mapping = get_seqmapping(args.inputfasta)
    # Invert the dictionary for rename_sequences
    mapping = {v:k for k,v in mapping.items()}
    # Rename all lines in each file using mapping
    for f in args.filestorename:
        rename_sequences(f, mapping)

def parse_args():
    parser = argparse.ArgumentParser(
       'Rename files from a seqview phyml run that did not get renamed back'
    )

    parser.add_argument(
        'inputfasta',
        help='Fast alignment file that was used to run phyml'
    )

    parser.add_argument(
        'filestorename',
        nargs='+',
        help='List of files to rename to original sequence names'
    )

    return parser.parse_args()
