from datetime import datetime
from config import root_dir, directories as dirs, rp_min_len
from common import ensure_dir, key_by_value, dump_buffer
from reporting import two_storey_bar_chart, run_FastQC
from seq_methods import merge_overlaps

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

def demux_illumina(dataset, max_pairs, run_id):
    """Demultiplex Illumina dataset.

    From separate forward/reverse read sets, combine read pairs and output
    to separate files for each sample based on barcode tags. As part of the
    process, reject read pairs that have mismatching tags or primers and trim
    the rest, removing primer+tag and low-quality sequences.
    """
    # identify inputs and outputs
    set_id = dataset['set_id']
    print " ", set_id
    run_root = root_dir+set_id+"/"+run_id+"/"
    ori_root = root_dir+set_id+"/"+dirs['master']
    fwd_file = ori_root+dataset['source_fwd']
    rev_file = ori_root+dataset['source_rev']
    demux_root = run_root+dirs['demux']
    report_root = run_root+dirs['reports']
    qc_dir = "qc_details/"
    qc_main_file = report_root+"quality_control.html"
    cntsplt = report_root+"sample_counts"
    ensure_dir(ori_root)
    ensure_dir(demux_root)
    ensure_dir(report_root)
    ensure_dir(report_root+qc_dir)
    # set up files for reporting
    html_comps = ["<p><b>Quality control for run "+run_id+"</b></p>",
                  "<p><img src='sample_counts.png' alt='sample_counts'/></p>",
                  "<p><table border='1'><tr>",
                  "<th>Sample</th>",
                  "<th>Accepted</th>",
                  "<th>Rejected</th>",
                  "<th>Total</th>",
                  "<th>% OK</th></tr>"]
    html_block = "".join(html_comps)
    open(qc_main_file, 'w').write(html_block)
    # prepare primers and barcodes info
    primers = dataset['primers']
    samples = dataset['samples']
    tag_pairs = samples.values()
    assert len(primers) >= 2
    assert len(samples) >= 1
    assert len(tag_pairs) >= 1
    # prepare container and files for output batching and reporting
    hits_dict = {}
    for sample_id in samples:
        hits_dict[sample_id] = {'buffer': [], 'countY': 0, 'countN': 0}
    # add containers for rejected read pairs
    hits_dict['bad_tags'] = {'buffer': [], 'countY': 0, 'countN': 0}
    hits_dict['bad_qual'] = {'buffer': [], 'countY': 0, 'countN': 0}
    # initialize files
    for sample_id in samples:
        dmx_out = demux_root+sample_id+"_readpairs.txt"
        open(dmx_out, 'w').write('')
    open(demux_root+"bad_tags"+"_readpairs.txt", 'w').write('')
    open(demux_root+"bad_qual"+"_readpairs.txt", 'w').write('')
    # iterate through reads
    pair_count = 0
    for titles, seqs, quals in FastqJointIterator(open(fwd_file),
                                                  open(rev_file)) :
        F_title = titles[0][0]
        R_title = titles[0][1]
        F_seq = seqs[0][0].upper()
        R_seq = seqs[0][1].upper()
        F_qual = quals[0][0]
        R_qual = quals[0][1]
        flip = False
        sample_id = False
        # iterate through barcode tags
        # TODO: implement more robust solution to ambiguous base problem
        for tag_pair in tag_pairs:
            L_tag1 = (tag_pair[0]+primers['fwdRA']).upper()
            L_tag2 = (tag_pair[0]+primers['fwdRG']).upper()
            R_tag = (tag_pair[1]+primers['rev']).upper()
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
                                else:
                                    R_clip = len(L_tag2)
                            else:
                                R_clip = len(L_tag1)
                            tag_hit = True
                            flip = True
                            F_clip = len(R_tag)
                            break
                else: # is there an L_tag in F_seq?
                    while True:
                        if not F_seq.find(L_tag1, 0, len(L_tag1)) is 0:
                            if not F_seq.find(L_tag2, 0, len(L_tag2)) is 0:
                                # no L_tag match -> reject
                                break
                            else:
                                F_clip = len(L_tag2)
                        else:
                            F_clip = len(L_tag1)
                        tag_hit = True
                        R_clip = len(R_tag)
                        break
                break
            if not tag_hit:     # continue iterating
                sample_id = False
            else:               # got it, stop iterating
                sample_id = key_by_value(samples, tag_pair)[0]
                break
        # in case no matches were found with any of the tags
        if not sample_id:
            sample_id = 'bad_tags'
        # for matched read pairs, clip off tag+primer and strip low qual runs
        else:
            F_trim = F_qual[F_clip:].find('##')
            if F_trim > -1:
                F_seq = F_seq[F_clip:F_clip+F_trim]
                F_qual = F_qual[F_clip:F_clip+F_trim]
            else:
                F_seq = F_seq[F_clip:]
                F_qual = F_qual[F_clip:]
            R_trim = R_qual[R_clip:].find('##')
            if R_trim > -1:
                R_seq = R_seq[R_clip:R_clip+R_trim]
                R_qual = R_qual[R_clip:R_clip+R_trim]
            else:
                R_seq = R_seq[R_clip:]
                R_qual = R_qual[R_clip:]
            if len(F_seq)+len(R_seq) < rp_min_len:
                # increment sample hit 'No' counter
                hits_dict[sample_id]['countN'] +=1
                sample_id = 'bad_qual'
        # bundle read data in ordered string
        readF = str("@%s\n%s\n+\n%s\n" % (F_title, F_seq, F_qual))
        readR = str("@%s\n%s\n+\n%s\n" % (R_title, R_seq, R_qual))
        if flip:
            read_pair = readR+readF
        else:
            read_pair = readF+readR
        # output to the appropriate buffer
        hits_dict[sample_id]['buffer'].append(read_pair)
        # increment sample 'Yes' hit counter
        hits_dict[sample_id]['countY'] +=1
        # when buffer capacity is reached, output to file and reset buffer
        if hits_dict[sample_id]['countY']% 100000==0:
            dmx_out = demux_root+sample_id+"_readpairs.txt"
            dump_buffer(dmx_out, hits_dict[sample_id]['buffer'])
            hits_dict[sample_id]['buffer'] = []
        # increment counter
        pair_count +=1
        # report on the progress
        if pair_count%1000000==0:
            print "\t", pair_count, "reads processed", datetime.now()
        if pair_count == max_pairs: # for testing purposes
            break
    print "\t", "Total", pair_count, "read pairs processed"
    print "\t", "Counts per sample:"
    # prepare graphing data containers
    pcntY = []
    pcntN = []
    sample_ids = []
    # write out whatever remains in each of the samples buffers
    for sample_id in samples:
        dmx_out = demux_root+sample_id+"_readpairs.txt"
        dump_buffer(dmx_out, hits_dict[sample_id]['buffer'])
        hits_dict[sample_id]['buffer'] = []
        acc = hits_dict[sample_id]['countY']
        rej = hits_dict[sample_id]['countN']
        print "\t\t", sample_id, acc, "pairs", datetime.now()
        pcntY.append(acc)
        pcntN.append(rej)
        sample_ids.append(sample_id)
        # generate FastQC report (use --noextract to not open zipped reports)
        run_FastQC(dmx_out, report_root+qc_dir, '--quiet', ' ')
        #print "see QC report"
        # add line in QC file
        link = qc_dir+sample_id+"_readpairs_fastqc/fastqc_report.html"
        html_comps = ["<tr>",
                      "<th><a href='"+link+"'>"+sample_id+"</a></th>",
                      "<td>", str(acc), "</td>",
                      "<td>", str(rej), "</td>",
                      "<td>", str(acc+rej), "</td>",
                      "<td>", str(int((float(acc)/(acc+rej))*100)),
                      "</td></tr>"]
        html_block = "".join(html_comps)
        open(qc_main_file, 'a').write(html_block)
    # write out whatever remains in the bad_qual buffer
    dmx_out = demux_root+"bad_qual_readpairs.txt"
    dump_buffer(dmx_out, hits_dict['bad_qual']['buffer'])
    hits_dict['bad_qual']['buffer'] = []
    print "\t\t", "rejected (low quality)", hits_dict['bad_qual']['countY'],\
    datetime.now()
    # generate FastQC report (use --noextract to not open zipped reports)
    run_FastQC(dmx_out, report_root+qc_dir, '--quiet', ' ')
    #print "see QC report"
    # add line in QC file
    link = qc_dir+"bad_qual_readpairs_fastqc/fastqc_report.html"
    html_comps = ["<tr>",
                  "<th><a href='"+link+"'>"+"bad_qual"+"</a></th>",
                  "<td>", '0', "</td>",
                  "<td>", str(hits_dict['bad_qual']['countY']), "</td>",
                  "<td>", str(hits_dict['bad_qual']['countY']), "</td>",
                  "<td>", '0',"</td></tr>"]
    html_block = "".join(html_comps)
    open(qc_main_file, 'a').write(html_block)
    # write out whatever remains in the bad_tags buffer
    dmx_out = demux_root+"bad_tags_readpairs.txt"
    dump_buffer(dmx_out, hits_dict['bad_tags']['buffer'])
    hits_dict['bad_tags']['buffer'] = []
    print "\t\t", "rejected (bad tags)", hits_dict['bad_tags']['countY'],\
    datetime.now()
    # generate FastQC report (use --noextract to not open zipped reports)
    run_FastQC(dmx_out, report_root+qc_dir, '--quiet', ' ')
    #print "see QC report"
    # add line in QC file
    link = qc_dir+"bad_tags_readpairs_fastqc/fastqc_report.html"
    html_comps = ["<tr>",
                  "<th><a href='"+link+"'>"+"bad_tags"+"</a></th>",
                  "<td>", '0', "</td>",
                  "<td>", str(hits_dict['bad_tags']['countY']), "</td>",
                  "<td>", str(hits_dict['bad_tags']['countY']), "</td>",
                  "<td>", '0',"</td></tr>"]
    html_block = "".join(html_comps)
    open(qc_main_file, 'a').write(html_block)
    # close table and add notes
    line_bq = "rejected after demultiplexing due to low sequence quality \
    (top stacks in bar chart)"
    line_bt = "could not be assigned to a sample due to mismatches in tag \
    and/or primer"
    html_comps = ["</table></p>",
                  "<p><b>", "bad_qual", ": </b>", line_bq,
                  "<br><b>", "bad_tags", ": </b>", line_bt, "</p>",]
    html_block = "".join(html_comps)
    open(qc_main_file, 'a').write(html_block)
    # add bad tags category for counts graphing (switch is on purpose)
    pcntY.append(hits_dict['bad_tags']['countN'])
    pcntN.append(hits_dict['bad_tags']['countY'])
    sample_ids.append('bad_tags')# check that the totals add up
    assert pair_count == sum(pcntY)+sum(pcntN)
    # plot the read counts per sample
    series = pcntY, pcntN
    legend = 'Accepted', 'Rejected'
    colors = 'g', 'r'
    titles = 'Number of read pairs', 'Read pairs per sample'
    two_storey_bar_chart(series, sample_ids, legend, colors, cntsplt, titles)

def merge_pair_libs(dataset, run_id):
    """Merge read pairs from Illumina sample libs and output FastA."""
    # identify inputs and outputs
    set_id = dataset['set_id']
    print " ", set_id
    run_root = root_dir+set_id+"/"+run_id+"/"
    dmx_root = run_root+dirs['demux']
    merged_root = run_root+dirs['merged']
    report_root = run_root+dirs['reports']
    master_file = run_root+dirs['merged']+run_id+".fas"
    ensure_dir(merged_root)
    ensure_dir(report_root)
    merger_file = report_root+"merged_pairs.html"
    cntsplt = report_root+"merge_counts"
    samples = dataset['samples']
    # set up files for reporting
    html_comps = ["<p><b>Read pairs merged for run ", run_id, "</b></p>",
                  "<p><img src='merge_counts.png' alt='merge_counts'/></p>",
                  "<p><table border='1'><tr>",
                  "<th>Sample</th>",
                  "<th>Accepted</th>",
                  "<th>Rejected</th>",
                  "<th>Total</th>",
                  "<th>% OK</th></tr>"]
    html_block = "".join(html_comps)
    open(merger_file, 'w').write(html_block)
    # initialize master file
    open(master_file, 'w').write('')
    # merge per sample (demuxed)
    merge_countA = []
    merge_countR = []
    sample_ids = samples.keys()
    for sample_id in sample_ids:
        print "\t", sample_id,
        lib_file = dmx_root+sample_id+"_readpairs.txt"
        merge_out = merged_root+sample_id+"_merged.fas"
        open(merge_out, 'w').write('')
        # prepare container and files for output batching and reporting
        buffer = []
        countY = 0
        countF = 0
        countN = 0
        # iterate through the read pairs
        count = 0
        for titles, seqs, quals in FastqGGIterator(open(lib_file)):
            count +=1
            seq1 = seqs[0]
            seq2 = seqs[1]
            qual1 = quals[0]
            qual2 = quals[1]
            # merge reads   TODO: better safeguard against merge failure
            try: merged = merge_overlaps(seq1, qual1, seq2, qual2)
            except: countF +=1
            else:
                if merged.find('N') > -1:
                    countN +=1  # if there are still N quality must be too low
                else:
                    countY +=1
                    # compose string for output
                    mcomps = [">",sample_id,"_",str(count),"\n",merged,"\n"]
                    mstring = "".join(mcomps)
                    # output to buffer
                    buffer.append(mstring)
            # when buffer capacity is reached, output to file and reset buffer
            if countY % 10000==0:
                dump_buffer(merge_out, buffer)
                dump_buffer(master_file, buffer)
                buffer = []
        # write out whatever remains in the buffer
        dump_buffer(merge_out, buffer)
        dump_buffer(master_file, buffer)
        # sum up
        assert countY+countF+countN == count
        print count, "pairs", datetime.now()
        print "\t\t", str(countY), "merged and accepted"
        print "\t\t", str(countN), "merged but rejected due to residual Ns"
        print "\t\t", str(countF), "failed to merge"
        # add line in QC file
        html_comps = ["<tr>",
                      "<th>", sample_id, "</b></th>",
                      "<td>", str(countY), "</td>",
                      "<td>", str(countN + countF), "</td>",
                      "<td>", str(count), "</td>",
                      "<td>", str(int((float(countY)/count)*100)),
                      "</td></tr>"]
        html_block = "".join(html_comps)
        open(merger_file, 'a').write(html_block)
        # pass values
        merge_countA.append(countY)
        merge_countR.append(countN+countF)
    # close table and add notes
    line_N = "either failed to merge or still contained Ns after merging"
    html_comps = ["</table></p>",
                  "<p><b>", "Rejected", ":</b> ", line_N, "</p>"]
    html_block = "".join(html_comps)
    open(merger_file, 'a').write(html_block)
    # plot the read counts per sample
    series = merge_countA, merge_countR
    lgnd = 'Accepted', 'Rejected'
    colors = 'g', 'r'
    titles = 'Number of read pairs', 'Read pairs merged per sample'
    two_storey_bar_chart(series, sample_ids, lgnd, colors, cntsplt, titles)
