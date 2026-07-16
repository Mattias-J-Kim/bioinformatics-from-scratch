# Bioinformatics From Scratch

Core bioinformatics algorithms implemented without BioPython,
to understand the underlying logic before relying on libraries.

## Why no BioPython?
(직접 채우기: 라이브러리에 의존하기 전에 알고리즘 자체를 이해하기 위함, 등)

## Modules
- **fasta_parser/** — multi-record FASTA 파일 파싱
- **gc_content/** — GC 비율 계산
- **reverse_complement/** — DNA 서열의 역상보서열 계산 (재귀)
- **orf_finder/** — 단일 frame에서 Open Reading Frame(ATG~stop codon) 탐색
- **six_frame_translation/** — 정방향/역방향 3-frame씩 총 6-frame 번역
- **hamming_distance/** — 서열 간 Hamming distance 계산 (⚠️ 현재 버그 있음, 파일 내 TODO 참고)

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
- [ ] hamming_distance — 버그 수정 필요 (전체 서열 순회 안 됨, index 시작점 오류)

## Notes
- `six_frame_translation/`은 자체 `reverse_complement()`를 내부에 갖고 있음 (반복문 버전, `reverse_complement/`의 재귀 버전과 별개 구현). 나중에 하나로 통합할지는 선택 사항.
