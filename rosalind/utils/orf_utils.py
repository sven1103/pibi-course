# Representation of all AA and their encoding triplets
__CODON_DICTIONARY = {
    'C': ['TGT', 'TGC'],
    'D': ['GAT', 'GAC'],
    'S': ['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT'],
    'Q': ['CAA', 'CAG'],
    'M': ['ATG'],
    'N': ['AAC', 'AAT'],
    'P': ['CCT', 'CCG', 'CCA', 'CCC'],
    'K': ['AAG', 'AAA'],
    'T': ['ACC', 'ACA', 'ACG', 'ACT'],
    'F': ['TTT', 'TTC'],
    'A': ['GCA', 'GCC', 'GCG', 'GCT'],
    'G': ['GGT', 'GGG', 'GGA', 'GGC'],
    'I': ['ATC', 'ATA', 'ATT'],
    'L': ['TTA', 'TTG', 'CTC', 'CTT', 'CTG', 'CTA'],
    'H': ['CAT', 'CAC'],
    'R': ['CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA'],
    'W': ['TGG'],
    'V': ['GTA', 'GTC', 'GTG', 'GTT'],
    'E': ['GAG', 'GAA'],
    'Y': ['TAT', 'TAC'],
    '*': ['TAG', 'TGA', 'TAA'],
}

# Complement bases
_COMPLEMENT = {
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C'
}


_RNA_ALPHABET = "AUCG"

_DNA_ALPHABET = "ATCG"


def valid_rna_letter(letter):
    return True if letter in _RNA_ALPHABET else False


def get_complement(s):
    return _COMPLEMENT[s]


def get_aminoacid(codon):
    for key in __CODON_DICTIONARY:
        if codon in __CODON_DICTIONARY[key]:
            return key
    return ""


def read_multi_fasta(fasta_file):
    all_sequences = []
    current_sequence = []
    with open(fasta_file, 'r') as file_handler:
        for line in file_handler.readlines():
            if '>' in line:
                if current_sequence:
                    all_sequences.append(*current_sequence)
                    current_sequence = []
                else:
                    continue
            else:
                current_sequence.append(line.rstrip())
    all_sequences.append(*current_sequence)
    return all_sequences
