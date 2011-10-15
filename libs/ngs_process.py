from config import directories as dirs
from Bio.SeqIO.QualityIO import FastqGeneralIterator

def FastqGGIterator(handle):
    """Iterate over 2 Fastq records as string tuples.

    Modified from Bio.SeqIO.QualityIO import FastqGeneralIterator to return
    two reads at a time.

    """
    handle_readline = handle.readline
    while True:
        line = handle_readline()
        if line == "" : return
        if line[0] == "@":
            break
    while True:
        title_lines = []
        seq_strings = []
        quality_strings = []
        count = 0
        while count < 2 :
            if line[0] != "@":
                raise ValueError("Bad formatting of record start lines.")
            title_line = line[1:].rstrip()
            title_lines.append(title_line)
            seq_string = handle_readline().rstrip()
            seq_strings.append(seq_string)
            while True:
                line = handle_readline()
                if not line:
                    raise ValueError("End of file without quality info.")
                if line[0] == "+":
                    second_title = line[1:].rstrip()
                    if second_title and second_title != title_line:
                        raise ValueError("Seq and qual captions differ.")
                    break
                seq_string += line.rstrip() #removes trailing newlines
            if " " in seq_string or "\t" in seq_string:
                raise ValueError("Whitespace not allowed in the sequence.")
            seq_len = len(seq_string)
            quality_string = handle_readline().rstrip()
            quality_strings.append(quality_string)
            while True:
                line = handle_readline()
                if not line : break #end of file
                if line[0] == "@":
                    if len(quality_string) >= seq_len:
                        break
                quality_string += line.rstrip()
            if seq_len != len(quality_string):
                raise ValueError("Lengths of seq and qual values differs "
                                 " for %s (%i and %i)." \
                                 % (title_line, seq_len, len(quality_string)))
            count +=1
        yield (title_lines, seq_strings, quality_strings)
        if not line : return #StopIteration at end of file
    assert False, "Should not reach this line"

def FastqJointIterator(handleA, handleB):
    """Iterate over Fastq records from two separate files as string tuples.

    Modified from Bio.SeqIO.QualityIO import FastqGeneralIterator to read
    from two separate files in parallel. Useful for shuffle-combining
    Illumina read sets.

    """
    handleA_readline = handleA.readline
    handleB_readline = handleB.readline
    while True:
        lineA = handleA_readline()
        lineB = handleB_readline()
        if lineA == "" and lineB == "": return
        if lineA[0] == "@" and lineB[0] == "@":
            break
    while True:
        title_lines = []
        seq_strings = []
        quality_strings = []
        count = 0
        while count < 2 :
            if lineA[0] != "@" and lineB[0] != "@":
                raise ValueError("Bad formatting of record start lines.")
            title_lineA = lineA[1:].rstrip()
            title_lineB = lineB[1:].rstrip()
            title_lines.append((title_lineA, title_lineB))
            seq_stringA = handleA_readline().rstrip()
            seq_stringB = handleB_readline().rstrip()
            seq_strings.append((seq_stringA, seq_stringB))
            while True:
                lineA = handleA_readline()
                lineB = handleB_readline()
                if not lineA:
                    raise ValueError("End of file without quality info.")
                if lineA[0] == "+" and lineB[0] == "+":
                    second_titleA = lineA[1:].rstrip()
                    second_titleB = lineB[1:].rstrip()
                    if second_titleA and second_titleA != title_lineA:
                        raise ValueError("Seq and qual captions differ (A).")
                    if second_titleB and second_titleB != title_lineB:
                        raise ValueError("Seq and qual captions differ (B).")
                    break
                seq_stringA += lineA.rstrip() #removes trailing newlines
                seq_stringB += lineB.rstrip() #removes trailing newlines
            if " " in seq_stringA or "\t" in seq_stringA:
                raise ValueError("Whitespace not allowed in sequence (A).")
            if " " in seq_stringB or "\t" in seq_stringB:
                raise ValueError("Whitespace not allowed in sequence (B).")
            seq_lenA = len(seq_stringA)
            seq_lenB = len(seq_stringB)
            quality_stringA = handleA_readline().rstrip()
            quality_stringB = handleB_readline().rstrip()
            quality_strings.append((quality_stringA, quality_stringB))
            while True:
                lineA = handleA_readline()
                if not lineA and not lineB: break #end of file
                if lineA[0] == "@" and lineB[0] == "@":
                    if len(quality_stringA) >= seq_lenA:
                        break
                    if len(quality_stringB) >= seq_lenB:
                        break
                quality_stringA += lineA.rstrip()
                quality_stringB += lineB.rstrip()
            if seq_lenA != len(quality_stringA):
                raise ValueError("Lengths of seq and qual values differs "
                                 " for %s (%i and %i) (A)." \
                                 % (title_lineA, seq_lenA,
                                    len(quality_stringA)))
            if seq_lenB != len(quality_stringB):
                raise ValueError("Lengths of seq and qual values differs "
                                 " for %s (%i and %i) (B)." \
                                 % (title_lineB, seq_lenB,
                                    len(quality_stringB)))
            count +=1
        yield (title_lines, seq_strings, quality_strings)
        if not lineA and not lineB: return #StopIteration at end of file
    assert False, "Should not reach this line"

def combine_illumina(dataset):
    """Combine forward and reverse read sets from Illumina.

    Consolidate both read sets in a single FastQ file, with each forward read
    immediately followed by the corresponding reverse read. This is meant to
    replace the "shuffle" Perl script provided with the Velvet assembler.

    """
    print " ", dataset['run_id']
    # identify inputs and outputs
    fwd_file = dataset['source_fwd']
    rev_file = dataset['source_rev']
    combined_file = dirs['combined_dir']+dataset['run_id']+".txt"
    # start shuffling
    combined_out = open(combined_file, 'w')
    read_count = 0
    for titles, seqs, quals in FastqGeneralIterator(open(fwd_file)) :
        print titles, seqs, quals
        read_count +=1
        if read_count >= 10:
            break