from __future__ import print_function
from rosalind.utils import orf_utils
import argparse

# Create the CLI
cli_parser = argparse.ArgumentParser(description=
                                     'Compute the Hamming distance matrix of '
                                     'given sequences.')
cli_parser.add_argument('multi_fasta', metavar='multi.fasta',
                        help='A FASTA file')
# Parse the command and export to namespace
args = cli_parser.parse_args()


def hamming_distance(string, other_string):
    hamming_dist = 0.0
    letter_counter = 0
    for letter in string:
        if letter not in other_string[letter_counter]:
            hamming_dist += 1.0
        letter_counter += 1
    return hamming_dist / len(string)


def calculate_distance_matrix(string, string_list):
    distances = []
    for other_string in string_list:
        distances.append(hamming_distance(string, other_string))
    return distances


def main():
    multi_fasta = vars(args)['multi_fasta']
    all_sequences = orf_utils.read_multi_fasta(multi_fasta)
    all_distances = []
    for sequence in all_sequences:
        all_distances.append(calculate_distance_matrix(sequence, all_sequences))

    for distance_row in all_distances:
        print("\t\t".join(map(str, distance_row)))


main()
