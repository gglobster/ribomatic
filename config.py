### Config file ###
# TODO: setup script to generate config file

## Data identification
datasets = [{'source_fwd': 'Novophage5_1_sequence.txt',
             'source_rev': 'Novophage5_2_sequence.txt',
             'set_id': 'Novophage5',
             'date': '9-22-11',
             'primers': {'fwdRA': 'caacgcgaAgaaccttacc', # workaround for R
                         'fwdRG': 'caacgcgaGgaaccttacc',
                         'rev': 'acaacacgagctgacgac'},
             'samples': {'CT34': ('tgca', 'actgc'),
                         'CT35': ('tgca', 'agcta'),
                         'CT36': ('tgca', 'acagtc'),
                         'CT37': ('tgca', 'tcgacg'),
                         'CT38': ('tgca', 'gcgac'),
                         'CT39': ('tgca', 'tacgt'),
                         'CT40': ('tgca', 'gtagtg'),
                         'CT42': ('tgca', 'gtca'),
                         'CT43': ('tgca', 'tact'),
                         'CT44': ('tgca', 'tcat'),
                         'CT45': ('tgca', 'tgca'),
                         'PM9': ('tgca', 'agt'),
                         'PM10': ('tgca', 'cga'),
                         'PM11': ('tgca', 'tac')}}]

## Directory structure

# project root directory
root_dir = 'data/'

directories = {
'master': 'master/',                # original data + master fasta file
'demux': 'demux/',                  # combined fwd/rev read sets
'merged': 'merged/',                # merged read pairs
'reports': 'reports/',              # reports
'otus': 'otus/',                    # OTU files
}

# Parameters

# read pair length min threshold (ensures overlap)
rp_min_len = 120