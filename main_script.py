from sys import argv, exit
from libs.common import ensure_dir
from libs.ngs_process import demux_illumina, merge_pair_libs
from libs.run_qiime import pick_otus, pick_rep_set, assign_taxonomy, \
    make_otu_table
from config import datasets, directories

print "\n", \
      "##################################################\n", \
      "### RiboMatic v. 0.1                           ###\n", \
      "### Copyright 2011 Geraldine A. Van der Auwera ###\n", \
      "##################################################\n", \

if len(argv) > 1 and argv[1] == '-h':
    print "Basic usage: \n", \
          "$ python main_script.py [step#]\n"
    exit()

if len(argv) < 2:
    step = 1
else:
    step = int(argv[1])

if step is 1:
    print "\n###", step, ". Demultiplex Illumina read sets ###"
    for dataset in datasets:
        demux_illumina(dataset)
    step +=1

if step is 2:
    print "\n###", step, ". Merge read pairs ###"
    for dataset in datasets:
        merge_pair_libs(dataset)
    step +=1

if step is 3:
    print "\n###", step, ". Cluster sequences into OTUs ###"
    for dataset in datasets:
        pick_otus(dataset)
    step +=1

if step is 4:
    print "\n###", step, ". Pick representative OTUs ###"
    for dataset in datasets:
        pick_rep_set(dataset)
    step +=1

if step is 5:
    print "\n###", step, ". Assign taxonomy to representative OTUs ###"
    for dataset in datasets:
        assign_taxonomy(dataset)
    step +=1

if step is 6:
    print "\n###", step, ". Make OTU table ###"
    for dataset in datasets:
        make_otu_table(dataset)
    step +=1

if step > 7:
    print "\n### Nothing more to do! ###\n"