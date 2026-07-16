def gc_content(seq):
    """
    서열의 GC 비율(%)을 계산. 소수점 둘째 자리까지 반올림.
    """
    gc_count = 0
    for base in seq:
        if base == "G" or base == "C":
            gc_count += 1
    gc_percent = gc_count * 100 / len(seq)
    return round(gc_percent, 2)


if __name__ == "__main__":
    print(gc_content("ATGCGCGTAGCTAGCTAGGCTATATGCGATCGAT"))
