import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import subprocess
from datetime import datetime
from config import root_dir, directories as dirs

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

def make_final_report(dataset):
    """Make final report HTML file."""
    # identify directories
    run_id = dataset['run_id']
    print " ", run_id
    run_root = run_id+"/"
    report_root = dirs['reports']
    heatmap_dir = report_root+"otu_heatmap/"
    comms_dir = report_root+"communities/"
    final_report_file = root_dir+run_root+run_id+"_report.html"
    # component links
    fqc_link = "http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/"
    qiime_link = "http://qiime.org"
    master_file = dirs['master']+run_id+".fas"
    qc_link = report_root+"quality_control.html"
    merge_link = report_root+"merged_pairs.html"
    stats_link = report_root+run_id+"_otu_table_stats.txt"
    heatmap_link = heatmap_dir+run_id+"_otu_table.html"
    area_link = comms_dir+"taxa_summary_plots/area_charts.html"
    bar_link = comms_dir+"taxa_summary_plots/bar_charts.html"
    # assemble report
    html_comps = ["<p><b>Final report for run "+run_id+"</b></p>",
                  "<p>Date compiled: "+str(datetime.now())+"</p>",
                  "<p>External components of the pipeline:<br><ul>",
                  "<li><a href='", fqc_link, "'>FastQC</a></li>",
                  "<li><a href='", qiime_link, "'>Qiime</a></li></ul></p>",
                  "<p>Master file of final 16S sequences:<br><ul>",
                  "<li><a href='", master_file, "'>"+run_id+".fas</a></li></ul></p>",
                  "<p>Report files:<br><ul>",
                  "<li><a href='", qc_link, "'>Quality control</a></li>",
                  "<li><a href='", merge_link, "'>Pair merging</a></li>",
                  "<li><a href='", stats_link, "'>Sample statistics</a></li>",
                  "<li><a href='", heatmap_link, "'>Heatmap of OTUs</a></li>",
                  "<li>Network of OTUs: see Qiime documentation</li>",
                  "<li><a href='", area_link,
                  "'>Area charts of taxonomic composition of samples</a></li>",
                  "<li><a href='", bar_link,
                  "'>Bar charts of taxonomic composition of samples</a></li>",
                  "</ul></p>"]
    html_block = "".join(html_comps)
    open(final_report_file, 'w').write(html_block)
    print "\t", "Final report complete at", run_root+run_id+"_report.html"