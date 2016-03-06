#!/usr/bin/python2.7
from __future__ import print_function
from rosalind.utils.errors import WrongAlphabetError
from rosalind.utils import orf_utils
import argparse
import sys

# Create the CLI
cli_parser = argparse.ArgumentParser(description=
                                     'Translate an RNA sequence into an amino'
                                     ' acid sequence.')
cli_parser.add_argument('rna_seq', metavar='UAGC....',
                        help='The RNA sequence string')
# Parse the command and export to namespace
args = cli_parser.parse_args()


def validate_rna(rna_seq):
    for letter in rna_seq:
        if not orf_utils.valid_rna_letter(letter):
            raise WrongAlphabetError(letter, "RNA")


def transform_to_dna(rna_seq):
    return rna_seq.replace("U", "T")


def main():
    rna_seq = vars(args)['rna_seq']
    try:
        validate_rna(rna_seq)
    except WrongAlphabetError as e:
        print("Wrong input sequence! {}".format(e.message))
        sys.exit(1)

    aa_seq = []
    for triplet in map(''.join, zip(*[iter(transform_to_dna(rna_seq))]*3)):
        aa_seq.append(orf_utils.get_aminoacid(triplet))

    print(''.join(aa_seq)[:-1])


main()

