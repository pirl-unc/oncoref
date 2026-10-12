"""Standard-code helpers for validating the historical nuclear ORF reference.

NCBI table 1: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1
Not a patient variant annotation engine.
"""

from itertools import product

CODONS = dict(
    zip(
        ("".join(c) for c in product("TCAG", repeat=3)),
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    )
)


def reverse_complement_dna(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def translate_reference(sequence):
    """Match the original catalogue builder's standard-code initiation rules."""
    residues = []
    for offset in range(0, len(sequence) - 2, 3):
        codon = sequence[offset : offset + 3]
        aa = "M" if offset == 0 and codon in ("ATG", "CTG", "TTG") else CODONS[codon]
        if aa == "*":
            return "".join(residues), True
        residues.append(aa)
    return "".join(residues), False
