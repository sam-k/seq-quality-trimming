# seq-quality-trimming
Trimming and merging reads from .ab1 Sanger sequencing files.

Part of the [human diet quantification project](http://el.ladlab.org:8080/research/), which aims to establish a model for quantifying dietary components from DNA metabarcoding data.

### Trimming sequences
- Identify a segment with the highest quality using Kadane's algorithm.
- Trim low-quality ends of the sequence.

### Merging sequences
- Align the forward and reverse PCR products.
- Merge the two sequences.

### Mapping sequences
- BLAST the merged sequence.