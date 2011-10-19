from sys import argv, exit
from libs.common import ensure_dir
from libs.ngs_process import demux_illumina, merge_pair_libs
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
    step = 0
else:
    step = int(argv[1])

if step is 0:
    ### STEP 0: Ensure that all base directories exist ###
    print "\n###", step, ". Set up the work environment ###"
    for dir_name in directories.keys():
        ensure_dir(directories[dir_name])
    step +=1

if step is 1:
    ### STEP 1: Demux Illumina read sets ###
    print "\n###", step, ". Demultiplex Illumina read sets ###"
    for dataset in datasets:
        demux_illumina(dataset)
    step +=1

if step is 2:
    ### STEP 2: Merge read pairs ###
    print "\n###", step, ". Merge read pairs ###"
    for dataset in datasets:
        merge_pair_libs(dataset)
    step +=1

if step > 3:
    print "\n### Nothing more to do! ###\n"