import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import subprocess
from common import ensure_dir
from config import root_dir, directories as dirs, rp_min_len
from os import path

def log_start_run(dataset, run_id, timestamp):
    """Record run initiated in the dataset log."""
    set_id = dataset['set_id']
    #print " ", set_id
    set_log = root_dir+set_id+"/"+set_id+"_log.html"
    param_file = run_id+"/"+dirs['reports']+run_id+"_parameters.txt"
    header = "<p><b>Log of processing runs for "+set_id+"</b></p><p><ul>"
    set_log_htm = ["<li><b>", run_id, "</b>&nbsp;- initiated ", timestamp,
                   " (<a href='", param_file, "'>parameters</a>)</li>"]
    linkline = "".join(set_log_htm)
    if not path.isfile(set_log): # first time running this dataset?
        open(set_log, 'w').write(header)
    open(set_log, 'a').write(linkline)

def log_end_run(dataset, run_id, timestamp):
    """Record run in the dataset log."""
    set_id = dataset['set_id']
    #print " ", set_id
    set_log = root_dir+set_id+"/"+set_id+"_log.html"
    run_report = run_id+"/"+run_id+"_report.html"
    set_log_htm = ["<li><b>", run_id, "</b>", "&nbsp;- completed ", timestamp,
                   " (<a href='", run_report, "'>report</a>)</li>"]
    linkline = "".join(set_log_htm)
    open(set_log, 'a').write(linkline)

def save_parameters(dataset, max_pairs, run_id, timestamp):
    """Save a copy of the dataset-specific parameters to file."""
    set_id = dataset['set_id']
    print " ", set_id
    run_root = root_dir+set_id+"/"+run_id+"/"
    report_root = run_root+dirs['reports']
    param_file = report_root+run_id+"_parameters.txt"
    ensure_dir(report_root)
    # primer data
    primers = dataset['primers']
    primers_list = []
    for primer_ID in primers:
        primers_list.append("\t".join([primer_ID, primers[primer_ID]]))
    primers_str = "\n".join(primers_list)
     # samples + barcode data
    samples = dataset['samples']
    samples_list = []
    for sample_ID in samples:
        samples_list.append("\t".join([sample_ID,
                                       samples[sample_ID][0],
                                       samples[sample_ID][1]]))
    samples_str = "\n".join(samples_list)
    # text block
    txt = ["# Run ID", run_id,
           "# Date generated", dataset['date'],
           "# Date processing initiated", timestamp,
           "# Special processing parameters",
           "# read pair length min threshold (ensures overlap)",
           str(rp_min_len),
           "# max number of read pairs to process",
           str(max_pairs),
           "# Dataset specifications",
           "# Illumina FastQ master files",
           dataset['source_fwd'], dataset['source_rev'],
           "# Amplification primers", primers_str,
           "# Sample ID\tLeft tag\tRight tag", samples_str]
    # write to file
    open(param_file, 'w').write("\n".join(txt))
    print "\t", "Run parameters saved to file"

def two_storey_bar_chart(series, labels, lgd_items, colors, fname, titles):
    """Draw bar chart of values with labels."""
    ind = np.arange(len(labels))    # the x locations for the groups
    width = 0.35                    # the width of the bars
    plt.subplot(111)
    plot1 = plt.bar(ind+width, series[0], width, color=colors[0])
    plot2 = plt.bar(ind+width, series[1], width, color=colors[1],
                    bottom=series[0])
    plt.xticks(ind+width, labels, rotation=30, size='small')
    plt.legend((plot1[0], plot2[0]), lgd_items)
    plt.ylabel(titles[0])
    plt.title(titles[1])
    plt.savefig(fname)
    plt.clf()

def run_FastQC(infile, outdir, quietness, unzip):
    """Make BLAST database from FASTA input file."""
    cline = "fastqc "+ infile +" -o "+ outdir +" "+ quietness+" "+unzip
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE, shell=True)
    output, error = child.communicate()
    return output

def make_final_report(dataset, run_id, start_stamp, end_stamp):
    """Make final report HTML file."""
    # identify directories
    set_id = dataset['set_id']
    print " ", set_id
    run_root = run_id+"/" # special case for local html linking
    report_root = dirs['reports']
    heatmap_dir = report_root+"otu_heatmap/"
    comms_dir = report_root+"communities/"
    final_report_file = root_dir+set_id+"/"+run_root+run_id+"_report.html"
    # component blurbs
    fastqc_blurb = "&nbsp;- A quality control tool for high throughput sequence data"
    qiime_blurb = "&nbsp;- Quantitative Insights Into Microbial Ecology"
    master_blurb = "Original dataset FastQ files"
    demux_blurb = "FastQ files of demultiplexed read pairs separated by sample"
    merged_blurb = "FastA files of accepted merged read pairs separated by sample + master file of all samples together"
    otus_blurb = "Text files and logs of OTU picking and processing by Qiime"
    reports_blurb = "Detailed reports for each pipeline step"
    qc_blurb = "- reject read pairs that have mismatches in the tag+primer region or have unacceptably low quality scores"
    merge_blurb = "- combine fwd and rev reads, taking best-scoring base where there is a mismatch"
    stats_blurb = "- relative amounts of accepted merged read pairs per sample"
    heatmap_blurb = "- relative abundance of OTUs per sample (interactive)"
    network_blurb = "- see Qiime documentation on how to use visualize these results"
    area_blurb = "- area charts (interactive)"
    bar_blurb = "- bar charts (interactive)"
    # component links
    runlog_link = dirs['reports']+run_id+"_parameters.txt"
    set_log_link = "../"+set_id+"_log.html"
    fqc_link = "http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/"
    qiime_link = "http://qiime.org"
    qc_link = report_root+"quality_control.html"
    merge_link = report_root+"merged_pairs.html"
    stats_link = report_root+run_id+"_otu_table_stats.txt"
    heatmap_link = heatmap_dir+run_id+"_otu_table.html"
    area_link = comms_dir+"taxa_summary_plots/area_charts.html"
    bar_link = comms_dir+"taxa_summary_plots/bar_charts.html"
    d_master = dirs['master']
    d_demux = dirs['demux']
    d_merged = dirs['merged']
    d_otus = dirs['otus']
    d_reports = dirs['reports']
    # code snippets
    li_a = "<li><a href='"
    e_li = "</li>"
    li_b = "<li><b>"
    b_nsp = "</b>&nbsp;"
    a_nsp = "</a>&nbsp;"
    # log files links
    run_log = "<a href='"+runlog_link+"'>run parameters log file</a>"
    set_log = "<a href='"+set_log_link+"'>Global log file</a>"
    # assemble report
    htm = ["<p><b>Run report for ", run_id, "</b></p>",
           "<p>Date/time initiated: ", start_stamp, "</p>",
           "<p>Date/time completed: ", end_stamp, "</p>",
           "<p>External components of the pipeline:<br><ul>",
           li_a, fqc_link, "'>FastQC</a>", fastqc_blurb, "</li>",
           "<li><a href='", qiime_link, "'>Qiime</a>",
           qiime_blurb, "</li></ul></p>",
           "<p>Directory structure of the results:<br><ul>",
           li_b, set_id,"/", d_master, b_nsp, master_blurb, e_li,
           li_b, set_id,"/", b_nsp, set_log, " for ", set_id, e_li,
           li_b, set_id,"/",run_id,"/", b_nsp, "This report & ",run_log, e_li,
           li_b, set_id,"/",run_id,"/", d_demux, b_nsp, demux_blurb, e_li,
           li_b, set_id,"/",run_id,"/", d_merged, b_nsp, merged_blurb, e_li,
           li_b, set_id,"/",run_id,"/", d_otus, b_nsp, otus_blurb, e_li,
           li_b, set_id,"/",run_id,"/", d_reports, b_nsp, reports_blurb, e_li,
           "</ul></p>",
           "<p>Pipeline steps / report files:<br><ul>",
           li_a, qc_link, "'>Demux + quality control", a_nsp, qc_blurb, e_li,
           li_a, merge_link, "'>Pair merging", a_nsp, merge_blurb, e_li,
           li_a, stats_link, "'>Sample statistics", a_nsp, stats_blurb, e_li,
           li_a, heatmap_link, "'>Heatmap of OTUs", a_nsp, heatmap_blurb, e_li,
           "<li>Network of OTUs", a_nsp, network_blurb, a_nsp, e_li,
           li_a, area_link, "'>Taxonomic composition of samples", a_nsp,
           area_blurb, e_li,
           li_a, bar_link, "'>Taxonomic composition of samples", a_nsp,
           bar_blurb, e_li,
           "</ul></p>"]
    html_block = "".join(htm)
    open(final_report_file, 'w').write(html_block)
    print "\t", "Run report complete"