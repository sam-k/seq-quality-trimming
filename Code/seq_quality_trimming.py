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
             "13E_004_P-12S-rev_BP006_A06",
             "15C_001_Abisporus_NS7_A06",
             "15C_002_Abisporus_ITS4_B06",
             "15C_003_Fvolutipes_NS7_C06",
             "15C_004_Fvolutipes_ITS4_D06",
             "15C_005_Htesselatus_NS7_E06",
             "15C_006_Htesselatus_ITS4_F06",
             "15C_007_Ledodes_NS7_G06",
             "15C_008_Ledodes_ITS4_H06",
             "15C_009_Peryngii_NS7_A07",
             "15C_010_Peryngii_ITS4_B07",
             "15C_011_Scerevisiae_NS7_C07",
             "15C_012_Scerevisiae_ITS4_D07",
             "15C_013_Rmicrosporus900_NS7_E07",
             "15C_014_Rmicrosporus900_ITS4_F07",
             "15C_015_Rmicrosporus1050_NS7_G07",
             "15C_016_Rmicrosporus1050_ITS4_H07"]

DNA_COMP = {"A":"T", "T":"A", "C":"G", "G":"C", "N":"N"}
# Find the reverse complement of a DNA sequence.
def rev_cmp(seq):
    revcmp_seq = ""
    for c in reversed(seq):
        revcmp_seq += DNA_COMP[c]
    return revcmp_seq

# All paired forward/reverse sequences.
# Each dict key is paired with a 2-tuple, corresponding to fwd/rev.
# str: [{names, reads, seqs, quals, trimmed_seqs}]
sequences = []
for _i in range(0, len(FILENAMES), 2):
    _names = (FILENAMES[_i], FILENAMES[_i+1])
    _reads = (SeqIO.read(WD + _names[0] + ".ab1", "abi"),
              SeqIO.read(WD + _names[1] + ".ab1", "abi"))
    _seqs = (str(_reads[0].seq), rev_cmp(str(_reads[1].seq)))
    _quals = ([x for x in _reads[0].letter_annotations["phred_quality"]],
              [x for x in _reads[1].letter_annotations["phred_quality"]][::-1])
    sequences.append({"names" : _names,
                      "reads" : _reads,
                      "seqs"  : _seqs,
                      "quals" : _quals})


    
""" trim by quality scores """

# Trim sequence by quality scores, using Kadane's algorithm.
# Cutoff is lowest acceptable quality score, used to normalize the scores around 0.
# Gap penalty is penalty per nt of including a low-quality base.
def trim_seq(seq, quality_scores, name=None, cutoff=None, gap_penalty=None):
    if cutoff is None:
        cutoff = np.average(quality_scores) * 0.9    # mean log, better cutoff metric needed
    if gap_penalty is None:
        gap_penalty = cutoff * -3 # arbitrary gap penalty, seems to work
    scores = [(x-cutoff if x > cutoff else x-cutoff+gap_penalty) for x in quality_scores]
    
    # Kadane's algorithm
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
        print("Entire sequence was low quality. Try a lower cutoff.")
        return ""
    
    max_ind = (ind_to_here if max_to_here > max_so_far else ind_so_far)
    trimmed_seq = seq[max_ind[0]:max_ind[1]+1]
    trimmed_scores = quality_scores[max_ind[0]:max_ind[1]+1]
    if name is not None:
        print(name)
    print("Trimmed seq: %d nt (%d-%d, %d%%)"
          % (len(trimmed_seq), max_ind[0]+1, max_ind[1]+1, 100*len(trimmed_seq)/len(seq)))
    print("Avg score: %.2f (std %.2f, rng %d-%d, >=%d, %d)\n"
          % (np.average(trimmed_scores), np.std(trimmed_scores),
             min(trimmed_scores), max(trimmed_scores), cutoff, gap_penalty))
    return trimmed_seq, max_ind

print("================== TRIMMING SEQUENCES ==================")

for _sq in sequences:
    _fwd_trim, _fwd_ind = trim_seq(_sq["seqs"][0], _sq["quals"][0], _sq["names"][0])
    _rev_trim, _rev_ind = trim_seq(_sq["seqs"][1], _sq["quals"][1], _sq["names"][1])
    _sq["trimmed_seqs"] = (_fwd_trim, _rev_trim)
    _sq["trim_indices"] = (_fwd_ind, _rev_ind)



""" merge forward/reverse sequences """

# Merge forward and reverse reads, using the Smith-Waterman algorithm.
# Different scores are assigned for matches, mismatches, and gaps in the alignment.
def merge_seqs(fwd, rev, name_fwd=None, name_rev=None, ind_fwd=None, ind_rev=None,
               match=None, mismatch=None, gap_penalty=None):
    # Build dictionary of scores for aligning bases
    dna_bases = ["A","T","C","G","N"]
    score_dict = {}
    for i in range(len(dna_bases)):
        for j in range(len(dna_bases)):
            score_dict[dna_bases[i]+dna_bases[j]] = match if i==j else mismatch
    score_dict["NN"] = mismatch
    
    # Initialize dynamic programming and traceback tables
    LEFT, DIAG, UP = range(3)   # define pointers: 0==Left, 1==Diagonal, 2==Up
    dp_table = [[0]*(len(rev)+1) for _ in range(len(fwd)+1)]    # base case = 0
    I = len(dp_table)
    J = len(dp_table[0])
    tb_table = [[[] for _ in range(len(rev)+1)] for _ in range(len(fwd)+1)] # no base case
    
    # Build tables
    max_align_score = 0
    max_back_index = (0, 0)
    for i in range(1, I):
        for j in range(1, J):
            scores = (gap_penalty + dp_table[i][j-1],  # left
                      score_dict[fwd[i-1]+rev[j-1]] + dp_table[i-1][j-1],   # diagonal
                      gap_penalty + dp_table[i-1][j],  # up
                      0)    # quit
            dp_table[i][j] = max_score = max(scores)
            if max_score != 0:
                for k in range(len(scores)):
                    if scores[k] == max_score:
                        tb_table[i][j].append(k)
                if max_score > max_align_score:
                    max_align_score = max_score
                    max_back_index = (i, j)
    
    # Traceback
    max_i, max_j = max_back_index
    max_front_index = (0, 0)
    fwd_align = rev_align = ""
    while len(tb_table[max_i][max_j]) > 0:
        max_front_index = (max_i, max_j)
        ptr = max(tb_table[max_i][max_j])
        if ptr != UP:
            if ptr == LEFT:
                fwd_align = "-" + fwd_align
            rev_align = rev[max_j-1] + rev_align
            max_j -= 1
        if ptr != LEFT:
            if ptr == UP:
                rev_align = "-" + rev_align
            fwd_align = fwd[max_i-1] + fwd_align
            max_i -= 1
    matches = 0
    indels = 0
    for i in range(len(fwd_align)):
        if fwd_align[i]==rev_align[i]:
            matches += 1
        elif fwd_align[i]=="-" or rev_align[i]=="-":
            indels += 1
    
    # Merge forward, aligned overlap, and reverse, and compile stats
    merged_seq = ""
    merged_nts = ""
    discarded_nts = [0, 0]
    if max_front_index[0] > max_front_index[1]:
        merged_seq += fwd[:max_front_index[0]]
        merged_nts += ("fwd %d-%d, "
                       % (ind_fwd[0], ind_fwd[0]+max_front_index[0]-1))
        discarded_nts[1] += max_front_index[1]-1
    else:
        merged_seq += rev[:max_front_index[1]]
        merged_nts += ("rev %d-%d, "
                       % (ind_rev[0], ind_rev[0]+max_front_index[1]-1))
        discarded_nts[0] += max_front_index[0]-1
    merged_seq += fwd[max_front_index[0]:max_back_index[0]]
    merged_nts += ("fwd %d-%d/rev %d-%d, "
                   % (ind_fwd[0]+max_front_index[0], ind_fwd[0]+max_back_index[0],
                      ind_rev[0]+max_front_index[1], ind_rev[0]+max_back_index[1]))
    if len(fwd)-max_back_index[0] < len(rev)-max_back_index[1]:
        merged_seq += rev[max_back_index[1]:]
        merged_nts += ("rev %d-%d"
                       % (ind_rev[0]+max_back_index[1]+1, ind_rev[1]+1))
        discarded_nts[0] += len(fwd)-max_back_index[0]
    else:
        merged_seq += fwd[max_back_index[0]:]
        merged_nts += ("fwd %d-%d"
                       % (ind_fwd[0]+max_back_index[0]+1, ind_fwd[1]+1))
        discarded_nts[1] += len(rev)-max_back_index[1]
    
    if name_fwd is not None and name_rev is not None:
        print(name_fwd + " + " + name_rev)
    overlap_len = max_back_index[0]-max_front_index[0]+1
    if overlap_len <= discarded_nts[0] or overlap_len <= discarded_nts[1]:
        print("Reads may not overlap or be too low quality. (%d nt, lost %d fwd/%d rev)\n"
              % (max_back_index[0]-max_front_index[0]+1, discarded_nts[0], discarded_nts[1]))
    else:
        print("Merged seq: %d nt (%s)" % (len(merged_seq), merged_nts))
        print("Overlap: %d nt (%.1f%% match, %d indels, lost %d fwd/%d rev)\n"
              % (max_back_index[0]-max_front_index[0]+1, 100*matches/len(fwd_align),
                 indels, discarded_nts[0], discarded_nts[1]))
    return merged_seq

print("================== MERGING SEQUENCES ===================")

for _sq in sequences:
    _sq["merged_seq"] = merge_seqs(_sq["trimmed_seqs"][0], _sq["trimmed_seqs"][1],
                                   _sq["names"][0], _sq["names"][1],
                                   _sq["trim_indices"][0], _sq["trim_indices"][1],
                                   2, -7, -7)   # arbitrary scores, seems to work