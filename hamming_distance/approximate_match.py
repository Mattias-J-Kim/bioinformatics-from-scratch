def approximate_match(text: str, pattern: str, d: int) -> list[int]:
    """
    Find all starting indices where `pattern` matches `text` with at most
    d mismatches (Hamming distance <= d).
 
    Unlike a simple GC-content sliding window, mismatch counts cannot be
    updated incrementally here: each shift changes which text position is
    compared against which pattern position, so every window is recomputed
    from scratch (O(n*k)).
 
    Args:
        text: sequence to search within.
        pattern: pattern to search for.
        d: maximum number of allowed mismatches.
 
    Returns:
        List of all starting indices with Hamming distance <= d.
    """
    k = len(pattern)
    matches = []
 
    for j in range(len(text) - k + 1):
        mismatches = 0
        for l in range(k):
            if text[j + l] != pattern[l]:
                mismatches += 1
        if mismatches <= d:
            matches.append(j)
 
    return matches
 
 
text = "GATATATGCATATACTT"
pattern = "ATAT"
d = 1
print(approximate_match(text, pattern, d))  # [1, 3, 9, 11]
