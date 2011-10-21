import time
from datetime import datetime
from sys import argv, exit
from libs.ngs_process import demux_illumina, merge_pair_libs
from libs.run_qiime import pick_otus, pick_rep_set, assign_taxonomy, \
    make_otu_table, make_qiime_reports
from libs.reporting import log_start_run, log_end_run, save_parameters, \
    make_final_report
from config import datasets

print "\n", \
      "##################################################\n", \
      "### RiboMatic v. 0.1                           ###\n", \
      "### Copyright 2011 Geraldine A. Van der Auwera ###\n", \
      "##################################################\n", \

if len(argv) > 1 and argv[1] == '-h':
    print "Basic usage: \n", \
          "$ python main_script.py [run_id] [step#] [max_pairs]\n", \
          "Note that these arguments are positional: order matters!\n"
    exit()

if len(argv) > 1:
    run_id = argv[1]
else:
    run_id = str(int(time.time())) # use timestamp as unique run identifier

if len(argv) > 2:
    max_pairs = int(argv[2])
else:
    max_pairs = 1000000000 # insanely large number = no limit

if len(argv) > 3:
    step = int(argv[3])
else:
    step = 0

start_timestamp = str(datetime.now())

if step is 0:
    print "\n###", step, ". Set up logging & reporting ###"
    for dataset in datasets:
        save_parameters(dataset, max_pairs, run_id, start_timestamp)
        log_start_run(dataset, run_id, start_timestamp)
    step +=1

if step is 1:
    print "\n###", step, ". Demultiplex Illumina read sets ###"
    for dataset in datasets:
        demux_illumina(dataset, max_pairs, run_id)
    step +=1

if step is 2:
    print "\n###", step, ". Merge read pairs ###"
    for dataset in datasets:
        merge_pair_libs(dataset, run_id)
    step +=1

if step is 3:
    print "\n###", step, ". Cluster sequences into OTUs ###"
    for dataset in datasets:
        pick_otus(dataset, run_id)
    step +=1

if step is 4:
    print "\n###", step, ". Pick representative OTUs ###"
    for dataset in datasets:
        pick_rep_set(dataset, run_id)
    step +=1

if step is 5:
    print "\n###", step, ". Assign taxonomy to representative OTUs ###"
    for dataset in datasets:
        assign_taxonomy(dataset, run_id)
    step +=1

if step is 6:
    print "\n###", step, ". Make OTU table ###"
    for dataset in datasets:
        make_otu_table(dataset, run_id)
    step +=1

if step is 7:
    print "\n###", step, ". Generate community analysis reports ###"
    for dataset in datasets:
        make_qiime_reports(dataset, run_id)
    step +=1

end_timestamp = str(datetime.now())

if step is 8:
    print "\n###", step, ". Compile final report and log run completion ###"
    for dataset in datasets:
        make_final_report(dataset, run_id, start_timestamp, end_timestamp)
        log_end_run(dataset, run_id, end_timestamp)
    step +=1

if step > 7:
    print "\n### Nothing more to do! ###\n"