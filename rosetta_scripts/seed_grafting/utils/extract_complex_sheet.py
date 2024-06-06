from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from utils import find_length
import argparse
import os
import pandas as pd
import numpy as np
import re
import sys
sys.path.append('..')
from global_vars import ch_target
from global_vars import ch_seed
from global_vars import dssp_path

################################################################################
# Utility script to extract and separate the target and sheet seed in two PDBs #
# Copyright A. Marchand (2024)                                                 #
################################################################################

def line_read_dssp(ding):
    """
    Takes a pdb filepath and returns name and DSSP secondary structure output
    """
    # Break up lines from file to get pdb and path
    pdb_filename = ding.split('/')[-1]
    pdb_filename = pdb_filename.replace("\n", "")
    pdb_id = pdb_filename[0:4]
    pdb_dir = ding.split('/')[:-1]
    pdb_dir = '/'.join([str(n) for n in pdb_dir])

    # Go to directory and load pdb to get DSSP
    os.chdir(pdb_dir)
    p = PDBParser()
    structure = p.get_structure(pdb_id, pdb_filename)
    model = structure[0]
    dssp = DSSP(model, pdb_filename, dssp=dssp_path)
    full_dssp = list(dssp)
    sec_str_result = [i[2] for i in full_dssp]
    sec_str_result = ''.join([str(n) for n in sec_str_result])

    return sec_str_result

def check_hetatm(file):
    with open(file,'r') as inf:
        datafile = inf.readlines()
    for line in datafile:
        if 'HETATM' in line:
            return True
    return False

base_dir = os.getcwd()
complex_file = sys.argv[1]
seed_id = sys.argv[2]

if check_hetatm(complex_file):
    het_adj = 1 #Adjustment for the presence of HETATM
else:
    het_adj = 0

dssp = line_read_dssp(complex_file)
os.chdir(base_dir)
start_seed, len_seed = find_length(complex_file, ch_seed)
dssp_seed = dssp[start_seed-1-het_adj:start_seed+len_seed-het_adj]
print(dssp_seed)

with open('./tmp/context_{}.pdb'.format(seed_id),'w') as file_context:
    with open(complex_file,'r') as in_pdb:
        for line in in_pdb:
            if (len(line.split())>3):
                if(((line[0:4] == "ATOM") or (line[0:6] == "HETATM")) and (line[21:22] == ch_target)):
                    file_context.write(line)
        file_context.write('TER')

previous_dssp = ''

with open('./tmp/seed_{}.pdb'.format(seed_id),'w') as file_seed:
    with open(complex_file,'r') as in_pdb:
        for line in in_pdb:
            if((line[0:4] == "ATOM") and (line[21:22] == ch_seed)):
                resid=int(line[22:26].strip(' '))
                if(dssp[resid-1-het_adj]=='E'):
                    file_seed.write(line)
                if(previous_dssp == 'E' and dssp[resid-1-het_adj] != 'E'):
                    file_seed.write('TER\n')
                previous_dssp = dssp[resid-1-het_adj]
