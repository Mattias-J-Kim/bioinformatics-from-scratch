# De Bruijn Graph Genome Assembler

A from-scratch implementation of short-read genome assembly: k-mers in, assembled
sequence out, with no graph libraries.

The general Eulerian path algorithm used here, together with the debugging log
and a second application of the same pattern in a non-biological domain, is
documented at
[`algorithms-practice/eulerian-path`](https://github.com/Mattias-J-Kim/algorithms-practice).
This README covers the biology.

---

## Reference

Compeau, P., Pevzner, P. & Tesler, G. (2011). "How to apply de Bruijn graphs to
genome assembly." *Nature Biotechnology* 29:987–991.

This is the structure underlying real short-read assemblers such as Velvet and
SPAdes.

---

## Construction

The counter-intuitive part of the De Bruijn formulation is that **k-mers are
edges, not nodes.** The nodes are (k−1)-mers, and each k-mer is the edge joining
its own prefix to its own suffix:

```
"ATG"  →  edge from "AT" to "TG"
```

```python
kmers = ["ATG", "TGC", "GCG"]

build_debruijn_graph(kmers)
# {"AT": ["TG"], "TG": ["GC"], "GC": ["CG"]}
```

Assembly is then the problem of walking every edge exactly once — an Eulerian
path — because every read must be used, and used only once.

The alternative formulation, treating reads as nodes and seeking a Hamiltonian
path, is the historically earlier approach and is NP-hard. The De Bruijn
reformulation is what makes assembly tractable at scale; this is the central
point of the Compeau et al. paper.

| Function | Responsibility |
|---|---|
| `build_debruijn_graph(kmers)` | k-mer list → adjacency list |
| `find_start_node(graph)` | node where `out_degree − in_degree == 1` |
| `eulerian_path(graph)` | Hierholzer traversal on a deep copy |
| `reconstruct_genome(path)` | overlap-collapse the node path into a string |
| `verify_assembly(genome, k, kmers)` | re-decompose and compare as a multiset |

---

## Results

**Test 1 — no repeated (k−1)-mers**

```
kmers    : ["ATG", "TGC", "GCG"]
assembled: ATGCG          ✓ matches source
```

**Test 2 — repeats present**

```
source   : TAATGCCATGGGATGTT
assembled: TAATGGGATGCCATGTT
identical to source : False
length correct      : True
verify_assembly     : True
```

**The assembled sequence differs from the source, and this is not a defect.**

The (k−1)-mer `AT` occurs at three separate positions in the source. In the
graph, all three collapse into a single node with three outgoing edges. Once that
happens, the information distinguishing "which `AT` was this" no longer exists —
it was destroyed when the genome was fragmented into k-mers, before the assembler
ever ran.

The graph therefore admits several distinct Eulerian paths, each of which
corresponds to a genome fully consistent with every observed k-mer. The traversal
returned one of them. No algorithm operating on this input can determine which
one was the original, because the input does not contain that information.

This is the repeat resolution problem, and it is the single largest obstacle in
real assembly. It is why the field moved toward paired-end reads (which constrain
the distance between two anchors) and long-read platforms such as PacBio and
Oxford Nanopore (which span repeat regions in a single read rather than
fragmenting across them).

**Consequence for testing.** Comparing the output to the source string is the
wrong correctness test, since it fails on valid assemblies. The implementation is
instead verified structurally:

- output length equals `len(kmers) + k − 1`
- re-decomposing the output into k-mers reproduces the input k-mer multiset

That is, the output is checked for consistency with the evidence actually
supplied, which is the only claim an assembler is entitled to make.

---

## Constraints Imposed

- No graph libraries (`networkx` excluded).
- Plain `dict` rather than `collections.defaultdict`, deliberately, in order to
  hit the adjacency-list initialization edge cases directly rather than have them
  handled invisibly.
- Start node determined by degree analysis rather than by picking an arbitrary
  key. Picking arbitrarily happens to work on Test 1 and fails on Test 2.
- The graph is deep-copied before traversal, since traversal consumes edges by
  popping them. A shallow copy still shares the inner destination lists, which
  would empty the caller's graph and make verification impossible.

---

## Open Items

- `find_start_node` raises `ValueError` when no node satisfies
  `out − in == 1`. That condition also holds for a fully balanced graph, which
  represents an Eulerian **circuit** — the correct structure for a **circular
  genome**. Bacterial chromosomes, plasmids and phage genomes such as PhiX174 are
  circular, so this implementation currently rejects a biologically common case
  it should handle. Fixing it means detecting the balanced case and choosing any
  node as the entry point.
- Real reads carry sequencing errors, which introduce spurious low-coverage
  branches ("tips" and "bubbles") into the graph. No error correction or graph
  cleaning is implemented; this assembler assumes perfect k-mers.
- Coverage is ignored entirely. Real assemblers use k-mer multiplicity to
  estimate copy number and to distinguish genuine repeats from errors.
- In-degree computation is O(V·E). A single pass over the edge list would reduce
  this.

---

## Note on Method

This problem was worked through with an AI assistant (Claude) used in a Socratic
mode: the assistant set the specification and test cases, and when I was stuck
responded with targeted questions and conceptual explanations rather than code.

All implementation code was written by me. Where the assistant identified a
defect, it described the failing condition and why it failed; the correction was
then made by me. The interpretation of the Test 2 result above reflects an
explanation I requested and then restated in my own terms.

This is stated explicitly because the value of this repository is the reasoning
trail, and a reasoning trail is only meaningful if its provenance is accurate.
