# Bioinformatics From Scratch

Core bioinformatics algorithms implemented without BioPython, to understand the
underlying logic before relying on libraries.

Python 3, standard library only. No external dependencies.

## Why no BioPython?

Implementing these algorithms manually — FASTA parsing, ORF detection, GC
content, translation, sequence comparison and assembly — before relying on
BioPython was a deliberate choice to understand the underlying logic rather than
treat it as a black box. This foundation supports later, more advanced work
(integrating statistical and optimization methods into sequence analysis
pipelines), and reflects the hands-on approach I'm building toward computational
genomics research.

## Modules

**Sequence handling**

- **fasta_parser/** — Multi-record FASTA file parsing
- **gc_content/** — GC content calculation, GC content sliding window (max GC window)
- **reverse_complement/** — Reverse complement of a DNA sequence (recursive)
- **orf_finder/** — Open Reading Frame (ATG–stop codon) detection in a single frame
- **six_frame_translation/** — Six-frame translation (3 forward + 3 reverse frames)

**Pattern and k-mer analysis**

- **hamming_distance/** — Hamming distance between sequences, approximate pattern matching
- **kmer_analysis/** — k-mer frequency counting, most frequent k-mer detection, clump finding (brute-force and optimized versions)
- **origin_of_replication/** — GC-skew based origin-of-replication (ori) candidate detection

**Alignment and assembly**

- **sequence_alignment/** — Global sequence alignment (Needleman–Wunsch style scoring: match/mismatch/gap), built on the Edit Distance DP recurrence
- **de_bruijn_assembly/** — De Bruijn graph construction and genome assembly via Eulerian path (Hierholzer's algorithm), with analysis of repeat-induced assembly ambiguity

### Pipelines

- **pipelines/** — Workflows combining the modules above into an end-to-end
  analysis. `pipeline2.py`: GC content profile across a genome using 100 bp
  non-overlapping windows (`fasta_parser` + `gc_content`).

### Data

- **data/** — Reference sequences used for testing.
  `NC_001422.1.fna`: bacteriophage phiX174 complete genome (5,386 bp),
  from NCBI RefSeq.

## Usage

Each module is a standalone script. Inputs are assigned as plain variables at the
bottom of the file rather than passed as arguments — edit those assignments and
run:

```
python <script_name>.py
```
Pipelines are run from the repository root:

```
python pipelines/pipeline2.py
```

## Related

The graph-algorithm reasoning behind `de_bruijn_assembly/`, together with a second
application of the same Eulerian path pattern in a non-biological domain, is
documented in
[algorithms-practice/eulerian-path](https://github.com/Mattias-J-Kim/algorithms-practice).

## Planned

- [ ] Circular genome support in `de_bruijn_assembly/` (Eulerian circuit case)
- [ ] Error-tolerant assembly: tip and bubble removal in the De Bruijn graph
- [ ] Local alignment (Smith–Waterman) alongside the existing global alignment
- [ ] Suffix array / BWT-based exact matching

## Notes

- **six_frame_translation/** has its own internal `reverse_complement()` (iterative
  version), implemented separately from the recursive version in
  `reverse_complement/`. Whether to merge them into one implementation later is
  optional.
- **sequence_alignment/** reuses the Edit Distance table shape; an early draft
  indexed both sequences with the same loop variable (`seq1[i-1]==seq2[i-1]`
  instead of `seq2[j-1]`), which didn't surface as a bug until tested with
  sequences of different lengths.
- **de_bruijn_assembly/** returns a sequence that differs from the source string
  when repeats are present. This is expected: repeated k-mers collapse into a
  single graph node, so multiple valid assemblies exist and the original cannot
  be recovered from k-mers alone. Correctness is therefore verified structurally
  (output length, and re-decomposition into the input k-mer multiset) rather than
  by string equality.
- **pipelines/pipeline2.py** discards the final partial window (86 bp of
  5,386). GC content from a shorter window is estimated from fewer bases and
  is not comparable to full-length windows on the same axis. phiX174 is
  circular, so wrapping the final window is the alternative — appropriate for
  GC-skew style analyses, unnecessary for a descriptive profile.

## Method

These modules were worked through with an AI assistant (Claude) used in a Socratic
mode: the assistant set specifications and test cases, and when I was stuck
responded with targeted questions and conceptual explanations rather than code.

All implementation code in this repository was written by me. Where the assistant
identified a defect, it described the failing condition and why it failed; the
correction was then made by me. 
One exception: the `__file__`-based path resolution in
`pipelines/pipeline2.py` was supplied by the assistant after I could not
construct the `os.path` calls unaided.

