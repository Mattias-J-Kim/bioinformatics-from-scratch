"""
Six-frame translation of a DNA sequence.

Translates a sequence in all three forward reading frames and all three
reverse-strand frames, using the standard genetic code.

Note:
    This module keeps its own iterative reverse_complement(), separate
    from the recursive version in reverse_complement/. They are not
    merged: importing across module folders would require path
    manipulation that the rest of the repository avoids outside of
    pipelines/.

    translate() stops at the first stop codon in each frame, so what is
    returned is the first peptide of each frame rather than the frame's
    full translation. A frame containing several ORFs will only report
    the first one.

    Codons that are not in the table — incomplete trailing codons, or
    codons containing ambiguity characters such as N — contribute nothing
    to the output. They are skipped silently rather than marked as X, so
    the returned peptide can be shorter than the frame implies.
"""

codon_table = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def reverse_complement(seq):
    """
    Return the reverse complement of a DNA sequence.

    Iterative implementation: each base is replaced by its complement,
    then the resulting string is reversed. Bases other than A, T, G and C
    are dropped, so the output can be shorter than the input.
    """
    comp = ""
    for base in seq:
        if base == "A":
            comp += "T"
        elif base == "T":
            comp += "A"
        elif base == "C":
            comp += "G"
        elif base == "G":
            comp += "C"
    return comp[::-1]


def translate(seq):
    """
    Translate a sequence in the three forward reading frames
    (offsets 0, 1, 2) and return the three peptides as a list.

    Each frame is read until its first stop codon, which is included in
    the output as "*".
    """
    results = []
    for offset in range(0, 3):
        protein = ""
        for i in range(offset, len(seq), 3):
            codon = seq[i:i + 3]
            aa = codon_table.get(codon, "")
            protein += aa
            if aa == "*":
                break
        results.append(protein)
    return results


def six_frame_translate(seq):
    """
    Translate all six reading frames: three on the given strand and three
    on its reverse complement.

    Returns:
        list of six peptides, forward frames first (offsets 0, 1, 2),
        then reverse-strand frames in the same offset order.
    """
    rc_seq = reverse_complement(seq)
    forward = translate(seq)
    reverse = translate(rc_seq)
    return forward + reverse


if __name__ == "__main__":
    print(six_frame_translate("AGTAGT"))
