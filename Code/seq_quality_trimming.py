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

# Biopython
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord



""" Initialize environment """

WD = os.path.expanduser("~/Projects/David Lab/Ind study projects/Seq trimming/Data/")
# Order by fwd/rev pairs. Even indices are fwd, odd indices are rev.
# Note that all rev sequences will be reverse-complemented to allow merging,
# which means all indices for rev reads will be mirrored.
FILENAMES = [
             "13E_001_BT-12S-fwd_BP005_F05",
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
             "15C_016_Rmicrosporus1050_ITS4_H07"
             ]
IMPORTED_FP = "seqs.fasta"
TRIMMED_FP  = "trimmed_seqs.fasta"
MERGED_FP   = "merged_seqs.fasta"
BLAST_FP    = "blast_results.txt"



# Import all paired forward/reverse sequences.
# Each dict key is paired with a 2-tuple, corresponding to fwd/rev.
# wd: Working directory, where reads are stored
# names: Filenames
# return: [{names, reads, seqs, quals, trimmed_seqs, trim_indices, merged_seq}]
def import_seqs(wd, filenames):
    sequences = []
    for i in range(0, len(filenames), 2):
        names = (filenames[i], filenames[i+1])
        reads = (SeqIO.read(wd + filenames[i] + ".ab1", "abi"),
                 SeqIO.read(wd + filenames[i+1] + ".ab1", "abi"))
        seqs  = (str(reads[0].seq), rev_cmp(str(reads[1].seq)))
        quals = ([x for x in reads[0].letter_annotations["phred_quality"]],
                  [x for x in reads[1].letter_annotations["phred_quality"]][::-1])
        sequences.append({"names":names, "seqs":seqs, "quals":quals})
    return sequences



# Find the reverse complement of a DNA sequence.
# seq: Sequence
# return: Rev cmp sequence
DNA_COMP = {"A":"T", "T":"A", "C":"G", "G":"C", "N":"N"}
def rev_cmp(seq):
    revcmp_seq = ""
    for c in reversed(seq):
        revcmp_seq += DNA_COMP[c]
    return revcmp_seq



# Trim sequence by quality scores, using Kadane's algorithm.
# seq: Raw sequence
# quality_scores: Phred quality scores for seq
# (Optional) name: Name of seq
# (Optional) fwd: Forward/reverse, used to position trimmed seq in original seq
# (Optional) cutoff: Lowest acceptable quality score, used to normalize scores around 0
# (Optional) low_penalty: Penalty per nt of including a low-quality base
# return: Trimmed seq, position in original seq
def trim_seq(seq, quality_scores, name=None, fwd=True, cutoff=None, low_penalty=None):
    if name is not None:
        print(name)
    if cutoff is None:
        cutoff = np.average(quality_scores) * 0.8    # mean log, better cutoff metric needed
    if low_penalty is None:
        low_penalty = cutoff * -3 # arbitrary gap penalty, seems to work
    scores = [(x-cutoff if x >= cutoff else x-cutoff+low_penalty) for x in quality_scores]
    
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
        print("Entire sequence was low quality. Try a lower cutoff.\n")
        return ("", "")
    max_ind = ind_to_here if max_to_here > max_so_far else ind_so_far
    trimmed_seq = seq[max_ind[0]:max_ind[1]+1]
    trimmed_scores = quality_scores[max_ind[0]:max_ind[1]+1]
    print("Trimmed seq: {} nt ({}-{}, {:.0f}%)".format(
            len(trimmed_seq),
            max_ind[0]+1 if fwd else len(seq)-max_ind[1],
            max_ind[1]+1 if fwd else len(seq)-max_ind[0],
            100*len(trimmed_seq)/len(seq)))
    print("Avg score: {:.2f} (std {:.2f}, rng {}-{}, >={:.0f}, {:.0f})\n".format(
            np.average(trimmed_scores), np.std(trimmed_scores),
            min(trimmed_scores), max(trimmed_scores), cutoff, low_penalty))
    return trimmed_seq, max_ind



# Merge forward and reverse reads, using the Smith-Waterman algorithm.
# Different scores are assigned for matches, mismatches, and gaps in the alignment.
# fwd/rev: Forward/reverse trimmed reads
# (Optional) name_fwd/rev: Names of fwd/rev reads
# (Optional) ind_fwd/rev: Positions of trimmed reads in original seqs
# (Optional) seq_fwd/rev: Original fwd/rev sequences
# (Optional) match/mismatch/gap_penalty: Scores for alignment
# return: merged seq
def merge_seqs(fwd, rev, name_fwd=None, name_rev=None, ind_fwd=None, ind_rev=None,
               len_seq_rev=None, match=5, mismatch=-5, gap_penalty=-5):
    if name_fwd is not None and name_rev is not None:
        print(name_fwd + " + " + name_rev)
    
    # Build dictionary of scores for aligning bases
    dna_bases = ["A","T","C","G","N"]
    score_dict = {}
    for i in range(len(dna_bases)):
        for j in range(len(dna_bases)):
            score_dict[dna_bases[i]+dna_bases[j]] = match if i==j else mismatch
    score_dict["NN"] = mismatch # no use matching unknown with unknown
    
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
    matches, indels = 0, 0
    for i in range(len(fwd_align)):
        if fwd_align[i]==rev_align[i]:
            matches += 1
        elif fwd_align[i]=="-" or rev_align[i]=="-":
            indels += 1
    
    # Merge forward, aligned overlap, and reverse, and compile stats
    merged_seq = ""
    merged_nts = ""
    discarded_nts = [0, 0]
    print_merge_stats = None not in (ind_fwd, ind_rev, len_seq_rev)
    if max_front_index[0] > max_front_index[1]: # before overlap
        merged_seq += fwd[:max_front_index[0]]
        if print_merge_stats:
            merged_nts += "fwd {}-{}, ".format(
                    ind_fwd[0]+1,
                    ind_fwd[0]+max_front_index[0]-1)
        discarded_nts[1] += max_front_index[1]-1
    else:
        merged_seq += rev[:max_front_index[1]]
        if print_merge_stats:
            merged_nts += "rev {}-{}, ".format(
                    len_seq_rev-ind_rev[0]+1,
                    len_seq_rev-ind_rev[0]-max_front_index[1])
        discarded_nts[0] += max_front_index[0]-1
    merged_seq += fwd[max_front_index[0]:max_back_index[0]] # overlap
    if print_merge_stats:
        merged_nts += "fwd {}-{}/rev {}-{}, ".format(
                ind_fwd[0]+max_front_index[0],
                ind_fwd[0]+max_back_index[0],
                len_seq_rev-ind_rev[0]-max_front_index[1]+1,
                len_seq_rev-ind_rev[0]-max_back_index[1]+1)
    if len(fwd)-max_back_index[0] < len(rev)-max_back_index[1]: # after overlap
        merged_seq += rev[max_back_index[1]:]
        if print_merge_stats:
            merged_nts += "rev {}-{}".format(
                    len_seq_rev-ind_rev[0]-max_back_index[1],
                    len_seq_rev-ind_rev[1])
        discarded_nts[0] += len(fwd)-max_back_index[0]
    else:
        merged_seq += fwd[max_back_index[0]:]
        if print_merge_stats:
            merged_nts += "fwd {}-{}".format(
                    ind_fwd[0]+max_back_index[0]+1,
                    ind_fwd[1]+1)
        discarded_nts[1] += len(rev)-max_back_index[1]
    
    overlap_len = max_back_index[0]-max_front_index[0]+1
    if overlap_len <= discarded_nts[0] or overlap_len <= discarded_nts[1]:
        merged_seq = None
        print("Reads may not overlap or be too low quality. ({} nt, lost {} fwd/{} rev)\n".format(
                max_back_index[0]-max_front_index[0]+1, discarded_nts[0], discarded_nts[1]))
    else:
        print("Merged seq: {} nt ({})".format(len(merged_seq), merged_nts))
        print("Overlap: {} nt ({:.1f}% match, {} indels, lost {} fwd/{} rev)\n".format(
                max_back_index[0]-max_front_index[0]+1, 100*matches/len(fwd_align),
                indels, discarded_nts[0], discarded_nts[1]))
    return merged_seq



# BLAST batch of reads.
# fp: FASTA file path
# return: BLAST records
def blast_seqs(fp):
    # Submit BLAST query
    return list(NCBIXML.parse(NCBIWWW.qblast("blastn", "nt", open(fp).read())))



# Write FASTA file.
# sequences: Dictionary of reads
# filename: Filename
# field: original, trimmed, or merged
def write_fasta(sequences, filename, field=""):
    records = []
    for sq in sequences:
        if field=="original":
            for k in range(2):
                records.append(SeqRecord(
                        Seq(sq["seqs"][k], IUPAC.IUPACAmbiguousDNA()),
                        id=sq["names"][k], description=""))
        elif field=="trimmed" or sq["merged_seq"] is None:
            for k in range(2):
                records.append(SeqRecord(
                        Seq(sq["trimmed_seqs"][k], IUPAC.IUPACAmbiguousDNA()),
                        id=sq["names"][k], description="({})".format(field)))
        elif field=="merged":
            records.append(SeqRecord(
                    Seq(sq["merged_seq"], IUPAC.IUPACAmbiguousDNA()),
                    id=", ".join(sq["names"]), description="(merged)"))
    SeqIO.write(records, filename, "fasta")



""" Main """

# Import fwd/rev reads
sequences = import_seqs(WD, FILENAMES)
write_fasta(sequences, IMPORTED_FP, "original")

# Trim fwd/rev reads
print("================== TRIMMING SEQUENCES ==================")
for _sq in sequences:
    _fwd_trim, _fwd_ind = trim_seq(_sq["seqs"][0], _sq["quals"][0], _sq["names"][0], True,
                                   20, -100)    # qual=20 == 99% accuracy
    _rev_trim, _rev_ind = trim_seq(_sq["seqs"][1], _sq["quals"][1], _sq["names"][1], False,
                                   20, -100)
    _sq["trimmed_seqs"] = (_fwd_trim, _rev_trim)
    _sq["trim_indices"] = (_fwd_ind,  _rev_ind)
write_fasta(sequences, TRIMMED_FP, "trimmed")

# Merge reads
print("================== MERGING SEQUENCES ===================")
for _sq in sequences:
    _sq["merged_seq"] = merge_seqs(_sq["trimmed_seqs"][0], _sq["trimmed_seqs"][1],
                                   _sq["names"][0],        _sq["names"][1],
                                   _sq["trim_indices"][0], _sq["trim_indices"][1],
                                   len(_sq["seqs"][1]),
                                   2, -7, -7)   # arbitrary scores, seems to work
write_fasta(sequences, MERGED_FP, "merged")

## BLAST merged reads and parse output
#print("================== BLASTING SEQUENCES ==================")
#blast_records = blast_seqs(MERGED_FP)
#_k = 0.0
#for _record in blast_records:
#    if _record.query.endswith("(unmerged)"):
#        if "blast" not in sequences[int(_k)]:
#            sequences[int(_k)]["blast"] = []
#            sequences[int(_k)]["blast_e"] = []
#        sequences[int(_k)]["blast"].append(_record.descriptions[0].title)
#        sequences[int(_k)]["blast_e"].append(_record.descriptions[0].e)
#        _k += 0.5
#    else:
#        sequences[int(_k)]["blast"] = _record.descriptions[0].title
#        sequences[int(_k)]["blast_e"] = _record.descriptions[0].e
#        _k += 1
#for _sq in sequences:
#    print("{} ({}merged)".format(
#            ", ".join(_sq["names"]), "un" if _sq["merged_seq"] is None else ""))
#    if type(_sq["blast"]) is list:
#        print("BLAST fwd: {} (e={:.1f})"  .format(_sq["blast"][0], _sq["blast_e"][0]))
#        print("BLAST rev: {} (e={:.1f})\n".format(_sq["blast"][1], _sq["blast_e"][1]))
#    else:
#        print("BLAST: {} (e={:.1f})\n".format(_sq["blast"], _sq["blast_e"]))