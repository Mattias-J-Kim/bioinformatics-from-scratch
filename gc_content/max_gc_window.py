def max_gc_window(seq: str, k: int) -> tuple[int, int]:
    """
    Find the window of length k with the highest GC content in a DNA sequence.
 
    Uses a sliding window: instead of recomputing the GC count for every
    window from scratch (O(n*k)), it updates the count incrementally by
    removing the leftmost base and adding the new rightmost base (O(n)).
 
    Args:
        seq: DNA sequence (A, T, G, C).
        k: window length.
 
    Returns:
        (start_index, gc_count) of the first window with the maximum GC count.
    """
    window = []
    gc_count = 0
    counts = []
 
    # initial window
    for i in range(k):
        window.append(seq[i])
        if seq[i] in ("C", "G"):
            gc_count += 1
    counts.append(gc_count)
 
    # slide window: drop leftmost base, add new rightmost base
    for j in range(len(seq) - k):
        if window[0] in ("C", "G"):
            gc_count -= 1
        window.pop(0)
        window.append(seq[k + j])
        if seq[k + j] in ("C", "G"):
            gc_count += 1
        counts.append(gc_count)
 
    best = max(counts)
    start_index = counts.index(best)
    return start_index, best
 
 
seq = "ATGCGCATGGCC"
k = 4
print(max_gc_window(seq, k))  # (2, 4)
