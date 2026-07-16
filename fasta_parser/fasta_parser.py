def parse_fasta(seq):
    """
    Multi-record FASTA 문자열을 파싱해서 {header: sequence} 딕셔너리로 반환.
    """
    lines = seq.split("\n")
    records = {}
    header = None
    for line in lines:
        if line.startswith(">"):
            header = line[1:]
            records[header] = ""
        else:
            records[header] += line
    return records


if __name__ == "__main__":
    example = """>seq1
ATGCGCGTAGCTAGCTA
GGCTATATGCGATCGAT
>seq2
ATATATATGCGCGCTTT
>seq3
GGGGCCCCAAAATTTT"""

    result = parse_fasta(example)
    for name, sequence in result.items():
        print(f"{name}: length={len(sequence)}")
