"""
Multi-record FASTA parsing.

Reads FASTA text into a {header: sequence} dictionary. Header lines are
stored without the leading ">", and a sequence split across several lines
is concatenated back into one string.
"""


def parse_fasta(text):
    """
    Parse multi-record FASTA text into a {header: sequence} dictionary.

    Args:
        text: contents of a FASTA file as a single string.

    Returns:
        dict mapping header (without ">") to the concatenated sequence.

    Raises:
        ValueError: if a sequence line appears before any header line, or
                    if the same header occurs twice.

    Note:
        Lines are stripped before use. Without this, a file saved with
        Windows line endings leaves a trailing "\\r" on every line, which
        ends up inside both the header and the sequence — the sequence is
        not rejected, just silently wrong, and its length is inflated by
        one per line.
    """
    records = {}
    header = None

    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            header = line[1:]
            if header in records:
                raise ValueError(f"Duplicate header: {header}")
            records[header] = ""
        else:
            if header is None:
                raise ValueError("Sequence line found before any header line")
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

    for name, sequence in parse_fasta(example).items():
        print(f"{name}: length={len(sequence)}")
