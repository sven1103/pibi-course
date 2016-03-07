#!/usr/bin/python2.7
from __future__ import print_function
from rosalind.utils.orf_utils import read_multi_fasta
import argparse

# Create the CLI
cli_parser = argparse.ArgumentParser(description=
                                     'Compute the Hamming distance matrix of '
                                     'given sequences.')
cli_parser.add_argument('multi_fasta', metavar='multi.fasta',
                        help='A FASTA file')
# Parse the command and export to namespace
args = cli_parser.parse_args()


# stores best result of shortest substring
# (the number of shortest substrings is optimally one = complete assembly)
best_substring = -1


def find_longest_common_substring(s1, s2):
    """ https://en.wikibooks.org/wiki/Algorithm_Implementation/
    Strings/Longest_common_substring
    :param s1: One read
    :param s2: other read
    :return: the longest common substring
    """
    m = [[0] * (1 + len(s2)) for i in range(1 + len(s1))]
    longest, x_longest = 0, 0
    for x in range(1, 1 + len(s1)):
        for y in range(1, 1 + len(s2)):
            if s1[x - 1] == s2[y - 1]:
                m[x][y] = m[x - 1][y - 1] + 1
                if m[x][y] > longest:
                    longest = m[x][y]
                    x_longest = x
            else:
                m[x][y] = 0
    return s1[x_longest - longest: x_longest]


def merge_reads(read, other_read, common_substring):
    """ Merges two reads at their common substring
    :param read: one read
    :param other_read: the other read
    :param common_substring: the common substring
    :return: the merged version of both reads
    """
    index_in_read = read.index(common_substring)
    index_in_other_read = other_read.index(common_substring)
    if index_in_read >= index_in_other_read:
        merged_reads = common_substring.join([read.replace(common_substring, ''),
                                             other_read.replace(common_substring, '')])
    else:
        merged_reads = common_substring.join([other_read.replace(common_substring, ''),
                                             read.replace(common_substring, '')])
    return merged_reads


def assemble_genome(read_list):
    """ Genome assembly
    :param read_list: The current read list (intermediate assembly status)
    :return: the assembled reads
    """
    # assembly complete, quit recursion
    if len(read_list) == 1:
        return read_list
    else:
        assembled_read_list = []
        already_merged = set([])
        for read in read_list:
            for other_read in read_list:
                if read not in already_merged and other_read not in already_merged:
                    longest_common_substring = find_longest_common_substring(read, other_read)
                    # merge reads, when overlap is greater half their size
                    if len(longest_common_substring) >= len(read)/2 and read is not other_read:
                        already_merged.add(read)
                        already_merged.add(other_read)
                        assembled_read_list.append(merge_reads(read, other_read, longest_common_substring))
                    else:
                        if read not in already_merged:
                            assembled_read_list.append(read)
                        if other_read not in already_merged:
                            assembled_read_list.append(other_read)

        assembled_read_list = filter(lambda _read: _read not in already_merged, assembled_read_list)

        if len(assembled_read_list) == len(read_list):
            # Assembly cannot be omptimized further, so return current assembly
            return assembled_read_list
        else:
            # Assembly was improved, continue optimization
            return assemble_genome(assembled_read_list)


def main():
    all_reads = read_multi_fasta(vars(args)['multi_fasta'])
    print("The best assembly from reads is:\n{0}".format(max(assemble_genome(all_reads))))

main()