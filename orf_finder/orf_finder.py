def find_orfs(sequence):
    """
    ATG(시작 코돈)부터 정지 코돈(TAA/TAG/TGA)까지의 ORF(Open Reading Frame)를
    3-nt 단위로 찾아서 리스트로 반환. (해당 frame 기준 단일 프레임 탐색)
    """
    orfs = []
    in_orf = False
    current = []
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i + 3]
        if not in_orf and codon == "ATG":
            in_orf = True
            current.append(codon)
        elif in_orf:
            current.append(codon)
            if codon in ("TAA", "TAG", "TGA"):
                orfs.append("".join(current))
                in_orf = False
                current = []
    return orfs


if __name__ == "__main__":
    # 실제 사용 시 fasta_parser.parse_fasta()로 만든 서열을 넣으면 됨
    # (모듈 간 import는 repo root 기준 실행 필요 — 여기선 독립 실행 예시만)
    example_seq = "ATGCGCGTAGCTAGCTAGGCTATATGCGATCGAT"
    orfs = find_orfs(example_seq)
    print(f"ORFs found={len(orfs)}")
    for orf in orfs:
        print(orf)
