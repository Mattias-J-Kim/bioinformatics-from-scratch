"""
Global Sequence Alignment (Needleman-Wunsch style scoring)
==============================================================

Learning note:
    Built directly on top of the Edit Distance recurrence from
    algorithms-practice/dynamic_programming — same table shape, same
    three-way choice (align / gap-in-seq1 / gap-in-seq2), just replacing
    "minimize edit operations" with "maximize alignment score".

    One indexing bug came up while implementing this independently:
        if seq1[i-1] == seq2[i-1]:   # wrong: reused i for both strings
    instead of:
        if seq1[i-1] == seq2[j-1]:   # correct: i indexes seq1, j indexes seq2
    This bug didn't throw an error or even produce a wrong score on the
    original test case, because seq1 and seq2 happened to have the same
    length (7), so i and j always advanced together in that specific
    example. It only surfaces as an IndexError once seq1 and seq2 have
    different lengths — a reminder that a matching accidental result on
    one test case doesn't confirm correct logic.

Problem:
    Given two DNA sequences, align them (allowing gaps to represent
    insertions/deletions) to maximize a scoring scheme:
        match    = +1
        mismatch = -1
        gap      = -2

    Example:
        seq1 = "GCATGCU"
        seq2 = "GATTACA"
        max score = -1

Relationship to Edit Distance:
    Edit Distance's three operations (substitute / insert / delete) map
    directly onto this problem's three choices (match-or-mismatch / gap
    in seq1 / gap in seq2). The only real difference is the objective:
    Edit Distance takes a MIN over operation counts; this problem takes
    a MAX over a scoring function.

DP table definition:
    dp[i][j] = best possible alignment score between seq1[0..i) and
               seq2[0..j)

Recurrence:
    dp[i][j] = max(
        dp[i-1][j-1] + (match_score if seq1[i-1]==seq2[j-1] else mismatch_score),
        dp[i-1][j] + gap_score,   # seq2 gets a gap
        dp[i][j-1] + gap_score,   # seq1 gets a gap
    )

Base cases:
    dp[i][0] = i * gap_score   (aligning the first i chars of seq1 against nothing)
    dp[0][j] = j * gap_score
"""


def sequence_alignment(seq1: str, seq2: str) -> int:
    MATCH, MISMATCH, GAP = 1, -1, -2
    n, m = len(seq1), len(seq2)
    dp = [[0] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        dp[i][0] = dp[i - 1][0] + GAP
    for j in range(1, m + 1):
        dp[0][j] = dp[0][j - 1] + GAP

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match_score = MATCH if seq1[i - 1] == seq2[j - 1] else MISMATCH
            dp[i][j] = max(
                dp[i - 1][j - 1] + match_score,
                dp[i - 1][j] + GAP,
                dp[i][j - 1] + GAP,
            )

    return dp[n][m]


seq1 = "GCATGCU"
seq2 = "GATTACA"
print(sequence_alignment(seq1, seq2))
