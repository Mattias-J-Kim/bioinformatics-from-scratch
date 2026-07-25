# Bioinformatics From Scratch

Core bioinformatics algorithms implemented without BioPython,
to understand the underlying logic before relying on libraries.

## Why no BioPython?
Implementing these algorithms manually — FASTA parsing, ORF detection, GC content,
translation, and sequence comparison — before relying on BioPython was a deliberate
choice to understand the underlying logic rather than treat it as a black box.
This foundation supports later, more advanced work (e.g. integrating statistical
methods like PCA and gradient-based optimization into sequence analysis pipelines),
and reflects the hands-on approach I'm building toward computational genomics research.

## Modules
- **fasta_parser/** — Multi-record FASTA file parsing
- **gc_content/** — GC content calculation, GC content sliding window (max GC window)
- **reverse_complement/** — Reverse complement of a DNA sequence (recursive)
- **orf_finder/** — Open Reading Frame (ATG–stop codon) detection in a single frame
- **six_frame_translation/** — Six-frame translation (3 forward + 3 reverse frames)
- **hamming_distance/** — Hamming distance calculation between sequences, approximate pattern matching
- **kmer_analysis/** — k-mer frequency counting, most frequent k-mer detection, clump finding (brute-force and optimized versions)
- **origin_of_replication/** — GC-skew based origin-of-replication (ori) candidate detection
- **sequence_alignment/** — Global sequence alignment (Needleman-Wunsch style scoring: match/mismatch/gap), built on the Edit Distance DP recurrence

## Usage
Inside each module folder:
```bash
python <script_name>.py <input_file>
```

## Status
- [x] fasta_parser
- [x] gc_content
- [x] reverse_complement
- [x] orf_finder
- [x] six_frame_translation
- [x] hamming_distance
- [x] kmer_analysis
- [x] origin_of_replication
- [x] sequence_alignment

## Notes
- `six_frame_translation/` has its own internal `reverse_complement()` (iterative version), implemented separately from the recursive version in `reverse_complement/`. Whether to merge them into one implementation later is optional.
- `sequence_alignment/` reuses the Edit Distance table shape; an early draft indexed both sequences with the same loop variable (`seq1[i-1]==seq2[i-1]` instead of `seq2[j-1]`), which didn't surface as a bug until tested with sequences of different lengths.
