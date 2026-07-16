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
- **fasta_parser/** — multi-record FASTA 파일 파싱
- **gc_content/** — GC 비율 계산
- **reverse_complement/** — DNA 서열의 역상보서열 계산 (재귀)
- **orf_finder/** — 단일 frame에서 Open Reading Frame(ATG~stop codon) 탐색
- **six_frame_translation/** — 정방향/역방향 3-frame씩 총 6-frame 번역
- **hamming_distance/** — 서열 간 Hamming distance 계산 

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

## Notes
- `six_frame_translation/`은 자체 `reverse_complement()`를 내부에 갖고 있음 (반복문 버전, `reverse_complement/`의 재귀 버전과 별개 구현). 나중에 하나로 통합할지는 선택 사항.
