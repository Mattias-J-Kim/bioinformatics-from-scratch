def reverse_complement(seq, index=None, _chars=None):
    """
    DNA 서열의 역상보서열(reverse complement)을 재귀로 계산.
    끝에서부터(index=len-1) 시작해서 앞으로 오면서 상보 염기로 치환.
    """
    if _chars is None:
        _chars = list(seq)
    if index is None:
        index = len(_chars) - 1
    if index < 0:
        return ""

    base = _chars[index]
    if base == "A":
        comp = "T"
    elif base == "T":
        comp = "A"
    elif base == "G":
        comp = "C"
    elif base == "C":
        comp = "G"
    else:
        raise ValueError(f"잘못된 문자: {base}")

    return comp + reverse_complement(seq, index - 1, _chars)


if __name__ == "__main__":
    n = input("염기서열을 입력해주세요: ")
    print(reverse_complement(n))
