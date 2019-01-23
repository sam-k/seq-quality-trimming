#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Trim seq by quality
DNA sequence
@date: January 17, 2019
@author: S. Kim
"""

import numpy as np
import os
from Bio import SeqIO  # Biopython's parser


""" set up environment """

wd = os.path.expanduser("~/Projects/David Lab/Ind study projects/Seq trimming/Data/")
filenames = ["13E_001_BT-12S-fwd_BP005_F05",
             "13E_002_BT-12S-rev_BP006_G05",
             "13E_003_P-12S-fwd_BP005_H05",
             "13E_004_P-12S-rev_BP006_A06"]

names = []  # str: [(forward_seq_file, reverse_seq_file)]
for _i in range(0, len(filenames), 2):
    names.append((filenames[_i], filenames[_i+1]))
seq_objects = [(SeqIO.read(wd + fwd + ".ab1", "abi"),
                SeqIO.read(wd + rev + ".ab1", "abi")) for (fwd, rev) in names]

# Paired forward/reverse sequences.
seqs = [(str(fwd.seq),
         str(rev.seq)) for (fwd, rev) in seq_objects]
# Quality values at each nucleotide. Parallel in structure with seqs list.
quals = [([x for x in fwd.letter_annotations["phred_quality"]],
          [x for x in rev.letter_annotations["phred_quality"]]) for (fwd, rev) in seq_objects]


""" trim by quality scores """

# Trim sequence by quality scores, using Kadane's algorithm.
# Threshold is lowest acceptable quality score, used to normalize the scores.
# Gap penalty is penalty per nt of including a low-quality base.
def trim_seq(seq, quality_scores, name=None, threshold=None, gap_penalty=None):
    if threshold is None:
        threshold = np.average(quality_scores) * 0.9
        gap_penalty = threshold * 3
    scores = [(x-threshold if x > threshold else x-threshold-gap_penalty) for x in quality_scores]
    
    max_to_here = 0
    max_so_far  = float("-inf")
    ind_to_here = [0, -1]
    ind_so_far  = [0, -1]
    for i in range(len(seq)):
        max_to_here += scores[i]
        ind_to_here[1] += 1
        if max_to_here > max_so_far:
            max_so_far = max_to_here
            ind_so_far = list(ind_to_here)
        if max_to_here < 0:
            max_to_here = 0
            ind_to_here = [i+1, i]
    if max_to_here <= 0 and max_so_far <= 0:
        print("Entire sequence was low quality. Try a lower threshold.")
        return ""
    
    max_ind = (ind_to_here if max_to_here > max_so_far else ind_so_far)
    trimmed_seq    = seq[max_ind[0]:max_ind[1]+1]
    trimmed_scores = quality_scores[max_ind[0]:max_ind[1]+1]
    if name is not None:
        print(name)
    print("Trimmed seq: %d nt (%d-%d, %d%%)"
          % (len(trimmed_seq), max_ind[0]+1, max_ind[1]+1, 100*len(trimmed_seq)/len(seq)))
    print("Avg score: %.2f (std %.2f, rng %d-%d, >=%d, -%d)\n"
          % (np.average(trimmed_scores), np.std(trimmed_scores),
             min(trimmed_scores), max(trimmed_scores), threshold, gap_penalty))
    return trimmed_seq

# Sequences with low-quality ends removed.
trimmed_seqs = []
for _i in range(len(seqs)):
    trimmed_seqs.append((trim_seq(seqs[_i][0], quals[_i][0], names[_i][0]),
                         trim_seq(seqs[_i][1], quals[_i][1], names[_i][1])))


""" merge forward/backward sequences """

