#!/usr/bin/python2.7
from __future__ import print_function

import argparse

from rosalind.utils import orf_utils

# Create the CLI
cli_parser = argparse.ArgumentParser(description=
                                     'Find all ORFs in a DNA sequence.')
cli_parser.add_argument('fasta_file', metavar='file.fasta',
                        help='A FASTA file containing the DNA sequence.')
# Parse the command and export to namespace
args = cli_parser.parse_args()


def _parse_fasta(file_name):
    sequence = []
    with open(file_name, 'r') as handler:
        for line in handler.readlines():
            if line.startswith('>'):
                continue
            else:
                sequence.append(line.rstrip())
    return ''.join(sequence)


def _get_orfs(sequence):
    return _parse_orf(sequence), _parse_orf(sequence[1:]),\
           _parse_orf(sequence[2:])


def _parse_orf(sequence):
    aa_sequence = []
    for triplet in map(''.join, zip(*[iter(sequence)]*3)):
        aa_sequence.append(orf_utils.get_aminoacid(triplet))
    return ''.join(aa_sequence)


def _get_reverse_complement(sequence):
    return ''.join(map(orf_utils.get_complement, sequence[::-1]))


def _evaluate_orf(sequence):
    if not sequence:
        return []
    # Split ORF at translation stop
    split_translation_stop = sequence.split("*") if sequence[:-1] is "*"\
        else sequence.split("*")[:-1]
    # Filter for split sequences that have a translation start signal ("M")
    filtered_list = []
    filtered_list.extend(filter(lambda x: "M" in x, split_translation_stop))

    # Now check how many methionine are in the filtered list strings
    # and make sub-orfs if necessary
    # i.e. MQRSMPENNN -> MQRSMPENNN and MPENNN
    final_list = []
    for orf in filtered_list:
        positions_met = [pos for pos, char in enumerate(orf) if char == "M"]
        map(lambda position: final_list.append(orf[position:]), positions_met)

    return final_list


def main():
    # Try to read the FASTA file
    sequence = _parse_fasta(vars(args)['fasta_file'])
    forward_orfs = list(_get_orfs(sequence))
    complement_orfs = list(_get_orfs(_get_reverse_complement(sequence)))
    result_list = []
    for seq in forward_orfs:
        result_list.extend(_evaluate_orf(seq))
    for seq in complement_orfs:
        result_list.extend(_evaluate_orf(seq))
    print(list(set(result_list)))


main()









