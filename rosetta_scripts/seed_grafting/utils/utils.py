import os,sys,string,gzip
import math
import Bio
from Bio.PDB import PDBParser

################################################################################
# A list of several utility scripts to prepare for the seed grafting protocol  #
# Copyright A. Marchand (2024)                                                 #
################################################################################

def find_chains(file):

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('struct', file)
    ch_list=[]
    for model in structure:
        for chain in model:
             ch_list.append(chain.get_id())
    return ch_list

def find_length( file, chain, frag_id = 1 ):

    list_resi=[]
    input_file = open(file)
    current_frag=1
    previous_line = []
    for line in input_file.readlines():
        split_line=line.split()
        if len(split_line) > 3 and split_line[0] == "ATOM" and split_line[4]==chain and current_frag==frag_id:
            list_resi.append(split_line[5])
        if line[0:3]=='TER' and previous_line[4] == chain :
            current_frag+=1
        previous_line = split_line

    length = int(list_resi[-1]) - int(list_resi[0]) + 1
    start = int(list_resi[0])

    return start,length

def find_min_dist(context_pdb, seed_pdb, resid, chain_target):

    dist=1000
    res_id=0

    input_pdb1 = open(seed_pdb.strip())

    for line in input_pdb1.readlines():
        if (len(line.split())>3):
            if((line[0:4] == "ATOM") and (line[21:22] != chain_target) and (line[22:26].strip(' ') == str(resid)) and (line[12:16].strip(' ') not in ['C','N','HA','CA','O','H'])):
                cs=[float(line[30:38].strip(' ')), float(line[38:46].strip(' ')), float(line[46:54].strip(' '))]
                input_pdb2 = open(context_pdb.strip())
                for line2 in input_pdb2.readlines():
                    if (len(line2.split())>3):
                        if(((line2[0:4] == "ATOM") or line2[0:6] == "HETATM")and (line2[21:22] == chain_target)):
                            ct=[float(line2[30:38].strip(' ')), float(line2[38:46].strip(' ')), float(line2[46:54].strip(' '))]
                            new_dist=math.sqrt(pow((cs[0]-ct[0]),2)+pow((cs[1]-ct[1]),2)+pow((cs[2]-ct[2]),2))

                            if new_dist<dist:
                                dist=new_dist
    return dist

def find_nfrag(seed_pdb):
    input_pdb = open(seed_pdb.strip())
    i=0
    for line in input_pdb.readlines():
        if(line[0:3]=='TER'):
            i+=1
    return i

def find_hotspots(context_pdb, seed_pdb, ch_seed, ch_target):

    outstr=''
    for j in range (1, find_nfrag(seed_pdb)+1,1):
        start, len_seed = find_length(seed_pdb, ch_seed, j)
        for i in range(start, start+len_seed, 1):
            dist = find_min_dist(context_pdb, seed_pdb, i, ch_target)
            if(dist < 2.5):
                outstr += str(i-start+1)
                outstr += ':'
        outstr = outstr.strip(':')
        outstr += ','

    outstr = outstr.strip(',')
    return outstr

def sanity_check(mol_id):
    error = False
    with open('./prep/{}.params'.format(mol_id)) as inf:
        for line in inf:
            split = line.split()
            if split[0] == 'ATOM':
                if line[4:6] != '  ':
                    print('WARNING: The following tabulation problem has been found in ./prep/{}.params file'.format(mol_id))
                    print(line.strip())
                    print('     ^')
                    error = True
            elif split[0] == 'BOND_TYPE':
                if line[9:11] != '  ':
                    print('WARNING: The following tabulation problem has been found in ./prep/{}.params file'.format(mol_id))
                    print(line.strip())
                    print('          ^')
                    error = True
            elif split[0] == 'ICOOR_INTERNAL':
                if line[14:18] != '    ':
                    print('WARNING: The following tabulation problem has been found in ./prep/{}.params file'.format(mol_id))
                    print(line.strip())
                    print('                 ^')
                    error = True
    if error:
        print('WARNING: Some formatting issues have been identified and can lead to unexpected behaviour when running RosettaScript.')
    else:
        print('Sanity check passed.')
