
"""
Hamming distance between two sequences of equal length.
"""


def hamming_distance(s1: str, s2: str) -> int:
    """
    Return the number of positions at which two equal-length strings
    s1 and s2 differ.

    Args:
        s1, s2: sequences of the same length.

    Returns:
        Number of differing positions.

    Raises:
        ValueError: if the two sequences differ in length. Hamming
                    distance compares position i against position i, so
                    it is undefined for sequences of different lengths —
                    the answer is not "the number of differences", it is
                    that the question doesn't apply.

    Example:
        hamming_distance("GGGCCGTTGGT", "GGACCGTTGAC") -> 3
    """
    if len(s1) != len(s2):
        raise ValueError(
            f"Sequences must be the same length: {len(s1)} != {len(s2)}"
        )

    distance = 0
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            distance += 1
    return distance


if __name__ == "__main__":
    print(hamming_distance("GGGCCGTTGGT", "GGACCGTTGAC"))  # 3
    print(hamming_distance("ATGC", "ATGC"))                # 0
