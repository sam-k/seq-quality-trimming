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

WD = os.path.expanduser("~/Projects/David Lab/Ind study projects/Seq trimming/")
ab1_fwd = SeqIO.read(WD+"Data/13E_001_BT-12S-fwd_BP005_F05.ab1", "abi")
ab1_rev = SeqIO.read(WD+"Data/13E_002_BT-12S-rev_BP006_G05.ab1", "abi")

seq_fwd = str(ab1_fwd.seq)
seq_rev = str(ab1_rev.seq)
qual_fwd = [x for x in ab1_fwd.letter_annotations["phred_quality"]]
qual_rev = [x for x in ab1_rev.letter_annotations["phred_quality"]]


""" trim by quality scores """

# Trim sequence by quality scores, using Kadane's algorithm.
# Threshold is lowest acceptable quality score.
# Gap penalty is penalty per nt of including a low-quality base.
def trim_seq(seq, quality_scores, threshold=0, gap_penalty=0):
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
    print("Trimmed seq: %d nt (%d-%d), Avg score: %.2f (std %.2f, rng %d-%d)"
          % (len(trimmed_seq), max_ind[0]+1, max_ind[1]+1,
             np.average(trimmed_scores), np.std(trimmed_scores), min(trimmed_scores), max(trimmed_scores)))
    return trimmed_seq

trimmed_seq_fwd = trim_seq(seq_fwd, qual_fwd, 40, 100)
trimmed_seq_rev = trim_seq(seq_rev, qual_rev, 40, 100)