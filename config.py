import numpy

### Config file ###
# TODO: setup script to generate config file

## Data identification
datasets = [{'source_fwd': 'Novophage5_1_sequence.txt',
             'source_rev': 'Novophage5_2_sequence.txt',
             'run_id': 'Novophage5',
             'date': '9-22-11'}]

## Directory structure

# project root directory
root_dir = 'test_data/'

directories = {
'ori_data_dir': root_dir+'original/',      # original data
'combined_dir': root_dir+'combined/',       # combined fwd/rev read sets

'blast_db_dir': root_dir+'blast_db/',      # blast databases
'reports_dir': root_dir+'reports/'         # reports
}