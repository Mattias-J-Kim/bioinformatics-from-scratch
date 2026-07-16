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
    주어진 서열을 3개 frame(offset 0,1,2)으로 번역해서 리스트로 반환.
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
    정방향 3-frame + 역상보 3-frame = 총 6개 frame 번역 결과를 반환.
    """
    rc_seq = reverse_complement(seq)
    forward = translate(seq)
    reverse = translate(rc_seq)
    return forward + reverse


if __name__ == "__main__":
    print(six_frame_translate("AGTAGT"))
