from datetime import datetime
from config import directories as dirs
from common import ensure_dir, key_by_value, dump_buffer
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
        if lineA == "" or lineB == "": return
        if lineA[0] == "@" or lineB[0] == "@": break
    while True:
        title_lines = []
        seq_strings = []
        quality_strings = []
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
                raise ValueError("End of file without quality info (A).")
            if not lineB:
                raise ValueError("End of file without quality info (B).")
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
            lineB = handleB_readline()
            if not lineA or not lineB: break #end of file
            if lineA[0] == "@":
                if len(quality_stringA) >= seq_lenA :
                    break
            quality_stringA += lineA.rstrip()
            if lineB[0] == "@":
                if len(quality_stringB) >= seq_lenB:
                    break
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
        yield (title_lines, seq_strings, quality_strings)
        if not lineA or not lineB: return #StopIteration at end of file
    assert False, "Should not reach this line"

def demux_illumina(dataset):
    """Demultiplex Illumina dataset.

    From separate forward/reverse read sets, combine read pairs and output
    to separate files for each sample based on barcode tags. As part of the
    process, reject read pairs that have mismatching tags or primers and trim
    the rest.

    """
    print " ", dataset['run_id']
    # identify inputs and outputs
    ori_root = dirs['ori_data']
    run_root = dataset['run_id']+"/"
    fwd_file = ori_root+run_root+dataset['source_fwd']
    rev_file = ori_root+run_root+dataset['source_rev']
    demux_root = dirs['demux']+run_root
    ensure_dir(demux_root)
    # prepare primers and barcodes info
    primers = dataset['primers']
    samples = dataset['samples']
    tags = samples.values()
    # prepare container for output batching and reporting
    hits_dict = {}
    for sample_id in samples:
        hits_dict[sample_id] = {'buffer': [], 'counter': 0}
    hits_dict['rejected'] = {'buffer': [], 'counter': 0}
    # iterate through reads
    pair_count = 0
    for titles, seqs, quals in FastqJointIterator(open(fwd_file),
                                                  open(rev_file)) :
        F_title = titles[0][0].lower()
        R_title = titles[0][1].lower()
        F_seq = seqs[0][0].lower()
        R_seq = seqs[0][1].lower()
        F_qual = quals[0][0].lower()
        R_qual = quals[0][1].lower()
        flip = False
        sample_id = False
        # iterate through barcode tags
        # TODO: implement more robust solution to ambiguous base problem
        for tag in tags:
            L_tag1 = (tag[0]+primers['fwdRA']).lower()
            L_tag2 = (tag[0]+primers['fwdRG']).lower()
            R_tag = (tag[1]+primers['rev']).lower()
            tag_hit = False
            while True:
                # start by checking For R_tag since there's only one
                if not R_seq.find(R_tag, 0, len(R_tag)) is 0:
                    if not F_seq.find(R_tag, 0, len(R_tag)) is 0:
                        # no R_tag match -> reject
                        break
                    else: # is there an L_tag in R_seq?
                        while True:
                            if not R_seq.find(L_tag1, 0, len(L_tag1)) is 0:
                                if not R_seq.find(L_tag2, 0, len(L_tag2)) is 0:
                                    # no L_tag match -> reject
                                    break
                            tag_hit = True
                            flip = True
                            break
                else: # is there an L_tag in F_seq?
                    while True:
                        if not F_seq.find(L_tag1, 0, len(L_tag1)) is 0:
                            if not F_seq.find(L_tag2, 0, len(L_tag2)) is 0:
                                # no L_tag match -> reject
                                break
                        tag_hit = True
                        break
                break
            if not tag_hit:     # continue iterating
                sample_id = False
            else:               # got it, stop iterating
                sample_id = key_by_value(samples, tag)[0]
                break
        # in case no matches were found with any of the tags
        if not sample_id:
            sample_id = 'rejected'
        # bundle read data in ordered string
        readF = str("@%s\n%s\n+\n%s\n" % (F_title, F_seq, F_qual))
        readR = str("@%s\n%s\n+\n%s\n" % (R_title, R_seq, R_qual))
        if flip:
            read_pair = readR+readF
        else:
            read_pair = readF+readR
        # output to the appropriate buffer
        hits_dict[sample_id]['buffer'].append(read_pair)
        # increment sample hit counter
        hits_dict[sample_id]['counter'] +=1
        # when buffer capacity is reached, output to file and reset buffer
        if hits_dict[sample_id]['counter']% 10000==0:
            dmx_out = demux_root+sample_id+"_readpairs.txt"
            dump_buffer(dmx_out, hits_dict[sample_id]['buffer'])
            hits_dict[sample_id]['buffer'] = []
            #print sample_id, "buffer reset @",
            # hits_dict[sample_id]['counter']
        # increment counter
        pair_count +=1
        # report on the progress
        if pair_count%100000==0:
            print "\t", pair_count, "reads processed", datetime.now()
#        if pair_count == 1000000: # for inspection purposes
#            break
    print "\t", "Counts per sample out of", pair_count, "total"
    # write out whatever remains in each of the samples buffers
    for sample_id in samples:
        dmx_out = demux_root+sample_id+"_readpairs.txt"
        dump_buffer(dmx_out, hits_dict[sample_id]['buffer'])
        hits_dict[sample_id]['buffer'] = []
        print "\t\t", sample_id, hits_dict[sample_id]['counter']
    # write out whatever remains in the rejected reads buffer
    dmx_out = demux_root+"rejected_readpairs.txt"
    dump_buffer(dmx_out, hits_dict['rejected']['buffer'])
    hits_dict['rejected']['buffer'] = []
    print "\t\t", "rejected", hits_dict['rejected']['counter']
