import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import subprocess

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
    # identify inputs and outputs
    