from Bio import PDB
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
import os,sys
sys.path.append('..')
from global_vars import dssp_path
import numpy as np
import os,sys
import pandas as pd

"""
This script calculates the number of atomic contact counts between the HELICAL binding seed and the small molecule
"""

def find_pdb_files():

    list_of_files=[]
    cmd = 'find ./out -name "*.pdb"'
    
    for file in os.popen(cmd).readlines():
        name = file[:-1]
        list_of_files.append(name)

    return list_of_files

def line_read_dssp(ding):
    """
    Takes a pdb filepath and returns name and DSSP secondary structure output
    """
    # Break up lines from file to get pdb and path
    pdb_filename = ding.split('/')[-1]
    pdb_filename = pdb_filename.replace("\n", "")
    pdb_dir = ding.split('/')[:-1]
    pdb_dir = '/'.join([str(n) for n in pdb_dir])
    cwd=os.getcwd()
    # Go to directory and load pdb to get DSSP
    os.chdir(pdb_dir)
    p = PDBParser()
    structure1 = p.get_structure(ding,pdb_filename)
    model = structure1[0]
    dssp = DSSP(model, pdb_filename, dssp=dssp_path)
    full_dssp = list(dssp)
    sec_str_result = [i[2] for i in full_dssp]
    sec_str_result = ''.join([str(n) for n in sec_str_result])
    os.chdir(cwd)
    return sec_str_result

def calculate_contacts(structure_file, chain_a_id, chain_b_id, drugname, SSE):
    structure = PDB.PDBParser().get_structure('protein', structure_file)
    dssp=line_read_dssp(structure_file)
    vdw_radius_lookup = {
    'H': 1.20,
    'C': 1.70,
    'N': 1.55,
    'O': 1.52,
    'Cl': 1.75,
    'P': 1.80,
    'S': 1.80
    }
    
    model = structure[0]
    chain_a = model[chain_a_id]
    chain_b = model[chain_b_id]

    contacts_drug = 0
    
    for atom_b in chain_b.get_atoms():
            
        for atom_a in chain_a.get_atoms():
            element_a = atom_a.element.strip() if atom_a.element else ''
            element_b = atom_b.element.strip() if atom_b.element else ''
            atom_parentid = atom_b.get_parent().id[1]
            
            atom_parent = atom_a.get_parent().resname
            if atom_parent!=drugname:
                continue

            vdw_radius_a = vdw_radius_lookup.get(element_a, 1.2)
            vdw_radius_b = vdw_radius_lookup.get(element_b, 1.2)

            threshold = vdw_radius_a + vdw_radius_b + 0.2
            distance = atom_a - atom_b

            if atom_parent == drugname and distance < threshold and dssp[atom_parentid-2] in SSE: 
                #-2 because DSSP string starts at index '0' and drug is not counted in DSSP
                contacts_drug += 1
                #print(atom_b.get_id(),atom_b.get_parent().id[1])
                break

    return contacts_drug

if len(sys.argv)<2:
    raise ValueError(f"Ligand identifier or SSE mode missing as an input")

mol=sys.argv[1]
sheet_only = sys.argv[2] == 'True'

if sheet_only == True:
    SSE=['E']
else:
    SSE=['H','B','E','G','I','T','S','-']

list_of_files=find_pdb_files()
dict_binders_prot={}
dict_binders_drug={}
i=1
for pdb in list_of_files:
    print("File %d out of %d done."%(i,len(list_of_files)))
    binder_name=os.path.basename(pdb).strip('.pdb')
    dict_binders_drug[binder_name]=calculate_contacts(pdb,'A','B',mol,SSE)
    i+=1

df_binders = pd.DataFrame.from_dict(dict_binders_drug, orient='index')
df_binders.index.name = 'design'
df_binders.rename(columns = {0:'mol_contacts'}, inplace = True)
df_binders.to_csv('contacts.csv')
print("Done.")
