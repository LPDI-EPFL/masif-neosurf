import os, sys
from utils import find_chains
sys.path.append('..')
from global_vars import ch_target
from global_vars import ch_seed

##################################################################################
# Utility script to extract and separate the target and helical seed in two PDBs #
# Copyright A. Marchand (2024)                                                   #
##################################################################################

complex_file = sys.argv[1]
seed_id = sys.argv[2]

for chain in find_chains(complex_file):
    if chain == ch_target:
        with open('./tmp/context_{}.pdb'.format(seed_id),'w') as file_context:
            with open(complex_file,'r') as in_pdb:
                for line in in_pdb:
                    if (len(line.split())>3):
                        if(((line[0:4] == "ATOM") or (line[0:6] == "HETATM")) and (line[21:22] == chain)):
                            file_context.write(line)
            file_context.write('TER')
        i=1
    elif chain == ch_seed:
        with open('./tmp/seed_{}.pdb'.format(seed_id),'w') as file_seed:
            with open(complex_file,'r') as in_pdb:
                for line in in_pdb:
                    if (len(line.split())>3):
                        if((line[0:4] == "ATOM") and (line[21:22] == chain)):
                            file_seed.write(line)
            file_seed.write('TER')
