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

WD = os.path.expanduser("~/Projects/David Lab/Ind study projects/Seq trimming/Data/")
FILENAMES = ["13E_001_BT-12S-fwd_BP005_F05",
             "13E_002_BT-12S-rev_BP006_G05",
             "13E_003_P-12S-fwd_BP005_H05",
             "13E_004_P-12S-rev_BP006_A06"]

# All paired forward/reverse sequences.
# Each dict key is paired with a 2-tuples, corresponding to fwd/rev.
# str: [{names, reads, seqs, quals, trimmed_seqs}]
sequences = []
for _i in range(0, len(FILENAMES), 2):
    names = (FILENAMES[_i], FILENAMES[_i+1])
    reads = (SeqIO.read(WD + names[0] + ".ab1", "abi"),
             SeqIO.read(WD + names[1] + ".ab1", "abi"))
    seqs = (str(reads[0].seq), str(reads[1].seq))
    quals = ([x for x in reads[0].letter_annotations["phred_quality"]],
             [x for x in reads[1].letter_annotations["phred_quality"]])
    sequences.append({"names" : names,
                      "reads" : reads,
                      "seqs"  : seqs,
                      "quals" : quals})


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

for _sq in sequences:
    _sq["trimmed_seqs"] = (trim_seq(_sq["seqs"][0], _sq["quals"][0], _sq["names"][0]),
                           trim_seq(_sq["seqs"][1], _sq["quals"][1], _sq["names"][1]))


""" merge forward/backward sequences """

