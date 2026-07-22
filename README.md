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

## Usage
각 모듈 폴더 안에서:
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

## Notes
- `six_frame_translation/`은 자체 `reverse_complement()`를 내부에 갖고 있음 (반복문 버전, `reverse_complement/`의 재귀 버전과 별개 구현). 나중에 하나로 통합할지는 선택 사항.
