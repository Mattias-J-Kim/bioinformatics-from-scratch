"""
Reverse complement of a DNA sequence.

Two implementations are kept here:

- reverse_complement_recursive : the recursion exercise. Python's default
  recursion limit (1000) means it fails with RecursionError on sequences
  longer than roughly 1,000 bp, so it cannot be used on
  data/NC_001422.1.fna (phiX174, 5,386 bp).
- reverse_complement : iterative version, used for real sequences. This is
  the implementation that previously lived inside six_frame_translation/;
  it was moved here to remove the duplication.
"""


def reverse_complement_recursive(seq, index=None, _chars=None):
    """
    Compute the reverse complement recursively.

    Starts at the last base (index = len - 1) and walks backwards,
    substituting each base for its complement.

    Note:
        Recursion depth grows with sequence length, so this fails on long
        sequences. Use reverse_complement() for anything genome-sized.
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
        raise ValueError(f"Invalid base: {base}")

    return comp + reverse_complement_recursive(seq, index - 1, _chars)


def reverse_complement(seq):
    """
    Compute the reverse complement iteratively. No length limit.

    Each base is replaced by its complement, then the resulting string is
    reversed.
    """
    comp = ""
    for base in seq:
        if base == "A":
            comp += "T"
        elif base == "T":
            comp += "A"
        elif base == "G":
            comp += "C"
        elif base == "C":
            comp += "G"
        else:
            raise ValueError(f"Invalid base: {base}")
    return comp[::-1]


if __name__ == "__main__":
    # Short sequence: check that both implementations agree.
    demo = "ATGCGCGTAGCTAGCTA"
    print("input :", demo)
    print("iter  :", reverse_complement(demo))
    print("rec   :", reverse_complement_recursive(demo))
    print("agree :", reverse_complement(demo) == reverse_complement_recursive(demo))
