def merge_overlaps(seq1, qual1, seq2, qual2):
    """Merge two sequences that overlap at the ends.

    This assumes both sequences are given in forward orientation. The second
    sequence is first complemented, then the two sequences are aligned to
    find the coordinates of overlap. The base calls are then compared for
    each position in the overlap; if there is a mismatch, compare the quality
    scores (reversed for seq2) using the ASCII values and keep the best one.

    >>> ord('B')
    66

    The alignment is done by iteratively scoring overlap configurations
    starting with a start-to-start case then pushing the second sequence
    further down along the first. Alignment is stopped and considered
    successful when the score is above 97% identity. If this point is not
    reached, uhhhhh. Not sure yet. We'll see.

    """
    # reverse-complement the second sequence
    seq2 = quick_revcom(seq2)
    # starting point and iterations depend on relative lengths
    if len(seq1) <= len(seq2):
        start = 0
        iter_max = len(seq1)
    else :
        start = len(seq1)-len(seq2)
        iter_max = len(seq2)
    iter_count = 0
    while iter_count < iter_max:
        print 'flebelebeleb'
        iter_count +=1

def quick_revcom(seq):
    """Take the DNA reverse complement of a sequence."""
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
    rc_seq = ''.join(rc_seq_temp)
    return rc_seq