import nwalign as nw

def merge_overlaps(seq1, qual1, seq2, qual2):
    """Merge two sequences that overlap at the ends.

    This assumes both sequences are given in forward orientation. The second
    sequence is first complemented, then the two sequences are aligned to
    find the coordinates of overlap. The base calls are then compared for
    each position in the overlap; if there is a mismatch, compare the quality
    scores (reversed for seq2) using the ASCII values and keep the best one.

    >>> ord('B')
    66

    The alignment is done using the nwalign package which implements the
    Needleman-Wunsch algorithm in C via Cython.

    """
    # reverse-complement the second sequence
    seq1 = seq1.upper()     # to make sure
    seq2 = quick_revcom(seq2)
    qual2 = qual2[::-1]
    # run the alignment
    alignment = nw.global_align(seq1, seq2, gap_open=-50, gap_extend=-2)
    seqln1 = alignment[0]
    seqln2 = alignment[1]
    # compose the merged sequence
    merge_list = []
    pos = 0
    qpos1 = 0
    qpos2 = 0
    while pos < len(alignment[0]):
        if seqln1[pos] is '-' or seqln1[pos] is 'N':
            merge_list.append(seqln2[pos])
            qpos1 +=1
        elif seqln2[pos] is '-' or seqln2[pos] is 'N':
            merge_list.append(seqln1[pos])
            qpos2 +=1
        elif seqln1[pos] is seqln2[pos]:
            merge_list.append(seqln1[pos])
        else:  # determine the consensus of the overlap using quality scores
            if ord(qual1[qpos1]) >= ord(qual2[qpos2]):
                merge_list.append(seqln1[pos])
            else:
                merge_list.append(seqln2[pos])
        pos +=1
        qpos1 +=1
        qpos2 +=1
    merged_seq = ''.join(merge_list)
    return merged_seq
        

def quick_revcom(seq):
    """Return the reverse complement of a DNA sequence.

    This is a low-overhead method for when we don't want to be dealing with
    the Bio package's feature-rich Seq objects.

    """
    rc_seq_temp = []
    for base in seq[::-1]:  # reverse
        if base is 'A':
            rc_seq_temp.append('T')
        elif base is 'T':
            rc_seq_temp.append('A')
        elif base is 'G':
            rc_seq_temp.append('C')
        elif base is 'C':
            rc_seq_temp.append('G')
        else:
            rc_seq_temp.append('N')
    rc_seq = ''.join(rc_seq_temp)
    return rc_seq