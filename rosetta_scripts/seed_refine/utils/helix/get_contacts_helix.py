from Bio import PDB
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

def calculate_contacts(structure_file, chain_a_id, chain_b_id, drugname):
    structure = PDB.PDBParser().get_structure('protein', structure_file)
    
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
            
            atom_parent = atom_a.get_parent().resname
            if atom_parent!=drugname:
                continue
            
            vdw_radius_a = vdw_radius_lookup.get(element_a, 1.2)
            vdw_radius_b = vdw_radius_lookup.get(element_b, 1.2)

            threshold = vdw_radius_a + vdw_radius_b + 0.2
            distance = atom_a - atom_b
            
            if atom_parent == drugname and distance < threshold:
                contacts_drug += 1
                #print(atom_b.get_id(),atom_b.get_parent().id[1])
                break

    return contacts_drug

if len(sys.argv)<1:
    raise ValueError(f"No ligand identifier was given as an input")

# Initialisation

mol=sys.argv[1]
list_of_files=find_pdb_files()
dict_binders_prot={}
dict_binders_drug={}
i=1

# Performing measurement of atomic contacts
for pdb in list_of_files:
    print("File %d out of %d done."%(i,len(list_of_files)))
    binder_name=os.path.basename(pdb).strip('.pdb')
    dict_binders_drug[binder_name]=calculate_contacts(pdb,'A','B',mol)
    i+=1

df_binders = pd.DataFrame.from_dict(dict_binders_drug, orient='index')
df_binders.index.name = 'design'
df_binders.rename(columns = {0:'mol_contacts'}, inplace = True)
df_binders.to_csv('contacts.csv')
print("Done.")
