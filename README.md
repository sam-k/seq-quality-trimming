# seq-quality-trimming
Trimming and merging reads from .ab1 Sanger sequencing files.

Part of the [human diet quantification project](http://el.ladlab.org:8080/), which aims to establish a model for quantifying dietary components from DNA metabarcoding data.

## Sequence processing
- Isolate a segment with the highest quality from the sequence.
- Trim low-quality ends of the sequence.
- Align and merge the forward/reverse PCR products.
- BLAST the merged sequence.