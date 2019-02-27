# seq-quality-trimming
Trimming and merging reads from .ab1 Sanger sequencing files.

Part of the [human diet quantification project](http://el.ladlab.org:8080/), which aims to establish a model for quantifying dietary components from DNA metabarcoding data.

## Quality trimming
- Isolate segment with the highest quality from the sequence, using Kadane's algorithm.
- Trim low-quality ends of the sequence.

## Merging reads
- Align the forward/reverse PCR products using the Smith-Waterman algorithm for local alignment.
- Merge the reads if the overlap is large enough.

## BLAST
- BLAST the merged sequence.