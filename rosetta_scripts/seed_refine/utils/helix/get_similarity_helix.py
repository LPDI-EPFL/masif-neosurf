import os
import glob
import gemmi
import datetime
import numpy as np
import pandas as pd
import Bio
from Bio import PDB
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.PDB import PDBParser

from scipy.spatial.distance import cdist

"""
This script calculates a similarity score for each helical seeds (0: Unique; 1: Identical with another seed)
The goal is to discard seeds that converge to the same hotspots residues
"""


def find_pdb_files():

    list_of_files=[]
    cmd = 'ls -d ./out/*.pdb'
    for file in os.popen(cmd).readlines():
        name = file[:-1]
        list_of_files.append(name)

    return list_of_files   

def get_hotspot_num(structure, chain_a_id, chain_b_id):

    van_der_waals_radii = {
    'H': 1.20,
    'C': 1.70,
    'N': 1.55,
    'O': 1.52,
    'Cl': 1.75,
    'P': 1.80,
    'S': 1.80
    }

    backbone_atoms = {'N', 'CA', 'C', 'O', 'H', 'HA'}

    contacts = []
    for residue_b in structure[0][chain_b_id]:
        for atom_b in residue_b:
            if residue_b.id[1] in contacts:
                break
            if atom_b.get_id() not in backbone_atoms:
                vdw_radius_b = van_der_waals_radii.get(atom_b.get_id(), 1.2)
                for residue_a in structure[0][chain_a_id]:
                    if residue_b.id[1] in contacts:
                        break
                    for atom_a in residue_a:
                        vdw_radius_a = van_der_waals_radii.get(atom_a.get_id(), 1.2)
                        combined_radius = vdw_radius_a + vdw_radius_b + 0.2
                        distance = atom_a - atom_b
                        if distance < combined_radius:
                            contacts.append(residue_b.id[1])
                            break
    return contacts

def get_hotspot_seq(structure, binder_chain,hotspot_num):
    hotspot_seq=''
    binder_residues=[]
    binder_seq=[]
    three_to_one = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }

    for residue in structure[0][binder_chain]:
        if PDB.is_aa(residue):
            binder_residues.append(residue.id[1])
            binder_seq.append(three_to_one.get(residue.get_resname()))

    for r,s in zip(binder_residues, binder_seq):
        if r in hotspot_num:
            hotspot_seq += s
        else:
            hotspot_seq += '-'

    return hotspot_seq


# Step 1: Initiatlization

list_of_files = find_pdb_files()
dict_binders = {}
binder_chain = 'B'
target_chain = 'A'
i = 0

#Step 2: Iteration on all PDBs

for pdb in list_of_files:
    i+=1
    parser = PDB.PDBParser()
    structure = parser.get_structure("protein", pdb)
    hotspot_num=get_hotspot_num(structure, target_chain, binder_chain)
    hotspot_seq = get_hotspot_seq(structure, binder_chain, hotspot_num)

    binder_name = os.path.basename(pdb).strip('.pdb')
    dict_binders[binder_name] = hotspot_seq
    print('File',i,'out of',len(list_of_files),'done.')

df_binders = pd.DataFrame.from_dict(dict_binders, orient='index')
df_binders.rename(columns = {0:'seq'}, inplace = True)
df_binders.index.name = 'design'

# Step 3: Score calculation and output

print("Performing alignement...")
for idx in df_binders.index:
    ali_score=0
    best_match=''
    for idx2 in df_binders.index:
        if idx!=idx2:
            seq1=df_binders.at[idx,'seq']
            seq2=df_binders.at[idx2,'seq']
            new_ali_score = pairwise2.align.globalms(seq1.strip("-"),seq2.strip("-"),1,0,-1,-1,score_only=True)
            if isinstance(new_ali_score, list) and len(new_ali_score) == 0:
                continue
            if new_ali_score>=ali_score:
                ali_score=new_ali_score
                best_match=str(idx2)
    seq=df_binders.at[idx,'seq'].strip('-')
    length=len(seq)
    if length == 0:
        uniqueness = 0
    else:
        uniqueness=ali_score/length #Score is normalize according to the sequence length
    df_binders.at[idx,'similarity']=uniqueness
    df_binders.at[idx,'best_match']=best_match

df_binders.to_csv('similarity.csv')
print("Done.")
