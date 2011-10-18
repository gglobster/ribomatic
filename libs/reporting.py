import numpy as np
import matplotlib
matplotlib.use('Agg')
#import matplotlib.pyplot as plt
import subprocess

def two_storey_bar_chart(series, labels, legend_items, colors, filename):
    """Draw bar chart of values with labels."""
    ind = np.arange(len(labels))    # the x locations for the groups
    width = 0.35                    # the width of the bars
#    plt.subplot(111)
#    plot1 = plt.bar(ind+width, series[0], width, color=colors[0])
#    plot2 = plt.bar(ind+width, series[1], width, color=colors[1],
#                    bottom=series[0])
#    plt.xticks(ind+width, labels)
#    plt.legend((plot1, plot2), legend_items)
#    plt.ylabel('Number of read pairs')
#    plt.title('Read pairs per sample')
#    plt.savefig(filename)

def run_FastQC(infile, outdir, quietness, unzip):
    """Make BLAST database from FASTA input file."""
    cline = "fastqc "+ infile +" -o "+ outdir +" "+ quietness+" "+unzip
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE, shell=True)
    output, error = child.communicate()
    return output