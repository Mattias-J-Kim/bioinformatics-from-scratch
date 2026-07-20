"""
origin_of_replication.py

Detects candidate origin-of-replication (ori) sites in a bacterial genome
using GC-skew, and finds DnaA-box-like binding motifs near those sites
using Hamming-distance-based approximate matching.

Background:
- Bacterial genomes are circular, and replication starts at ori and
  proceeds bidirectionally.
- The leading and lagging strands undergo C->T mutations at different
  rates during replication, which causes the ratio of G to C to flip
  around ori. This produces a measurable pattern called GC-skew.
- The position where skew reaches its minimum is a strong candidate
  for the ori location.
- Binding motifs such as DnaA boxes accumulate mutations over
  evolutionary time, so exact string matching is not sufficient.
  Approximate matching (allowing up to d mismatches) is needed to
  recover biologically meaningful occurrences.

Related: builds on the same "local pattern detection" theme as
find_clumps (see kmer_analysis/ in bioinformatics-from-scratch) -
that function detects clumps of exact-match k-mers, while this file
extends the idea to approximate (mismatch-tolerant) matching.
"""


def skew_array(genome: str) -> list:
    """
    Return the running (G count - C count) at every prefix of genome.

    The returned list has length len(genome) + 1, since it starts with
    the value 0 before any base has been read.

    Example:
        skew_array("CATG") -> [0, -1, -1, -1, 0]
    """
    skew = [0]
    current = 0
    for base in genome:
        if base == "G":
            current += 1
        elif base == "C":
            current -= 1
        # A and T leave the running value unchanged
        skew.append(current)
    return skew


def min_skew_positions(genome: str) -> list:
    """
    Return all indices where skew_array reaches its minimum value.
    These positions are candidate origin-of-replication (ori) sites.
    """
    skew = skew_array(genome)
    min_val = min(skew)  # compute once, reuse for every comparison
    return [i for i, val in enumerate(skew) if val == min_val]


def hamming_distance(s1: str, s2: str) -> int:
    """
    Return the number of positions at which two equal-length strings
    s1 and s2 differ.

    Example:
        hamming_distance("GGGCCGTTGGT", "GGACCGTTGAC") -> 3
    """
    distance = 0
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            distance += 1
    return distance


def approximate_pattern_matching(pattern: str, genome: str, d: int) -> list:
    """
    Return the starting positions of every substring of genome whose
    Hamming distance to pattern is at most d.

    A match does not need to be exact; up to d mismatches are allowed.
    This models the fact that real binding sites (e.g. DnaA boxes)
    accumulate mutations and rarely appear as perfect repeats.

    Example:
        approximate_pattern_matching("GA", "ATGACGA", 1) -> [2, 5]
    """
    k = len(pattern)
    positions = []
    for i in range(0, len(genome) - k + 1):
        window = genome[i:i + k]
        if hamming_distance(pattern, window) <= d:
            positions.append(i)
    return positions


if __name__ == "__main__":
    # skew_array check
    print(skew_array("CATG"))  # [0, -1, -1, -1, 0]

    # min_skew_positions check
    print(min_skew_positions("CATGGGCATCGGCCATACGCC"))  # [21]

    # hamming_distance check
    print(hamming_distance("GGGCCGTTGGT", "GGACCGTTGAC"))  # 3

    # approximate_pattern_matching check
    print(approximate_pattern_matching("GA", "ATGACGA", 1))  # [2, 5]

    # A slightly larger example, for comparison with the exact-match
    # clump-finding approach used in kmer_analysis/find_clumps.
    seq = (
        "CGGACTCGACAGATGTGAAGAACGAAAATGTGACGACGGACGACGGCATGAACGAAAAC"
        "TGACAGAAGTGACAAAGTGA"
    )
    print("min skew position:", min_skew_positions(seq))
