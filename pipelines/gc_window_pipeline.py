"""
GC content sliding-window pipeline.
 
This is the first script in the repository that connects previously
standalone modules into a single workflow:
 
    fasta_parser.parse_fasta  ->  {header: sequence} from a FASTA file
    gc_content.gc_content     ->  GC percentage of a sequence
 
For each record in the input file it reports the whole-sequence GC content,
then a GC profile computed over consecutive non-overlapping windows.
 
Test data
---------
Bacteriophage phiX174 complete genome
NCBI RefSeq accession NC_001422.1 (5,386 bp)
 
Design notes
------------
1. Path resolution
   File paths are built from __file__ rather than passed as relative paths.
   A relative path is resolved against the current working directory, which
   is where the command was typed, not where the script lives. Anchoring on
   __file__ makes the script behave identically regardless of the directory
   it is invoked from.
 
2. Handling the trailing partial window
   5,386 bp does not divide evenly into 100 bp windows; the final segment is
   only 86 bp. That remainder is dropped rather than reported. GC content is
   a ratio, so a short window still yields a value, but it is estimated from
   fewer bases and therefore varies more. Reporting it alongside full windows
   would place non-comparable values on the same axis. Discarding incomplete
   trailing windows is the common default in window-based genome analysis.
 
   The alternative would be to treat the genome as circular and wrap the
   final window around to the start, which is biologically correct for
   phiX174. That is worth doing for analyses where the origin of replication
   matters (e.g. GC skew), but it adds coordinate-reporting complexity that a
   descriptive GC profile does not need.
"""
import sys
import os
 
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(BASE_DIR, "fasta_parser"))
sys.path.append(os.path.join(BASE_DIR, "gc_content"))
 
from fasta_parser import parse_fasta
from gc_content import gc_content
 
DATA_PATH = os.path.join(BASE_DIR, "data", "NC_001422.1.fna")
 
with open(DATA_PATH) as f:
    content = f.read()
    a = parse_fasta(content)
 
WINDOW = 100
 
for header, seq in a.items():
    records = gc_content(seq)
    print(f"{header}")
    print(f"Overall GC: {records}%")
    # Incomplete trailing window is dropped: GC from a shorter window is
    # noisier and not comparable to full-length windows.
    for i in range(0, len(seq) - WINDOW + 1, WINDOW):
        print(f"{i}bp - {i+WINDOW}bp: {gc_content(seq[i:i+WINDOW])}%")
