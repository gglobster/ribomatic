import subprocess
from config import root_dir, directories as dirs
from common import ensure_dir

# TODO: adapt cline to run with the regular CLI qiime binary (here=macqiime)

def pick_otus(dataset):
    """Pick OTUs using Uclust with Qiime."""
    # identify inputs and outputs
    run_id = dataset['run_id']
    print " ", run_id
    run_root = root_dir+run_id+"/"
    otus_dir = run_root+dirs['otus']
    ensure_dir(otus_dir)
    master_file = run_root+dirs['master']+run_id+".fas"
    # run the command
    comps = ["macqiime", "pick_otus.py", "-i", master_file, "-o", otus_dir]
    cline = " ".join(comps)
    try:
        child = subprocess.Popen(str(cline), stdout=subprocess.PIPE,
                                 shell=True)
        output, error = child.communicate()
    except: raise
    else: print "\t", "OTUs picked"

def pick_rep_set(dataset):
    """Pick representative set of OTUs with Qiime."""
    # identify inputs and outputs
    run_id = dataset['run_id']
    print " ", run_id
    run_root = root_dir+run_id+"/"
    otus_dir = run_root+dirs['otus']
    master_file = run_root+dirs['master']+run_id+".fas"
    otus_file = otus_dir+run_id+"_otus.txt"
    rep_otus_file = otus_dir+run_id+"_rep_set.fas"
    # run the command
    comps = ["macqiime", "pick_rep_set.py", "-i", otus_file,
             "-f", master_file, "-o", rep_otus_file]
    cline = " ".join(comps)
    try:
        child = subprocess.Popen(str(cline), stdout=subprocess.PIPE,
                                 shell=True)
        output, error = child.communicate()
    except: raise
    else: print "\t", "OTUs representative set selected"

def assign_taxonomy(dataset):
    """Assign taxonomy to representative OTUs using RDP with Qiime."""
    # identify inputs and outputs
    run_id = dataset['run_id']
    print " ", run_id
    run_root = root_dir+run_id+"/"
    otus_dir = run_root+dirs['otus']
    rep_otus_file = otus_dir+run_id+"_rep_set.fas"
    # run the command
    comps = ["macqiime", "assign_taxonomy.py", "-i", rep_otus_file,
             "-o", otus_dir, "-m", "rdp"]
    cline = " ".join(comps)
    try:
        child = subprocess.Popen(str(cline), stdout=subprocess.PIPE,
                                 shell=True)
        output, error = child.communicate()
    except: raise
    else: print "\t", "Taxonomy assigned to representative OTUs"

def make_otu_table(dataset):
    """Make table of OTUs and summarize sample statistics with Qiime."""
    # identify inputs and outputs
    run_id = dataset['run_id']
    print " ", run_id
    run_root = root_dir+run_id+"/"
    otus_dir = run_root+dirs['otus']
    otus_file = otus_dir+run_id+"_otus.txt"
    tax_ass_file = otus_dir+run_id+"_rep_set_tax_assignments.txt"
    table_file = otus_dir+run_id+"_otu_table.txt"
    stats_file = otus_dir+run_id+"_otu_table_stats.txt"
    # run the command
    comps = ["macqiime", "make_otu_table.py", "-i", otus_file,
             "-t", tax_ass_file, "-o", table_file]
    cline = " ".join(comps)
    try:
        child = subprocess.Popen(str(cline), stdout=subprocess.PIPE,
                                 shell=True)
        output, error = child.communicate()
    except: raise
    else: print "\t", "OTU table compiled",
    # now summarize the stats
    comps = ["macqiime", "per_library_stats.py", "-i", table_file]
    cline = " ".join(comps)
    try:
        child = subprocess.Popen(str(cline), stdout=subprocess.PIPE,
                                 shell=True)
        output, error = child.communicate()
    except: raise
    else:
        open(stats_file, 'w').write(output)
        print "and stats summarized"
