from Bio import PDB
import numpy as np
import os,sys
import pandas as pd


################################################################################
# Utility script that returns all target residues in contact with a drug       #
# Copyright A. Marchand (2024)                                                 #
################################################################################

def contact_res(structure_file, chain_a_id, drugname):
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
    backbone_atoms = {'N', 'CA', 'C', 'O', 'H', 'HA'}
    contacts_drug = []

    for atom_a1 in chain_a.get_atoms():
        atom_parent1 = atom_a1.get_parent().resname
        if atom_parent1!=drugname:
            continue   
        for atom_a2 in chain_a.get_atoms():
            element_a1 = atom_a1.element.strip() if atom_a1.element else ''
            element_a2 = atom_a2.element.strip() if atom_a2.element else ''
            
            atom_parent2 = atom_a2.get_parent().resname
            if atom_parent2 == drugname or atom_a2.get_id() in backbone_atoms:
                continue
            
            vdw_radius_a1 = vdw_radius_lookup.get(element_a1, 1.2)
            vdw_radius_a2 = vdw_radius_lookup.get(element_a2, 1.2)

            threshold = vdw_radius_a1 + vdw_radius_a2 + 0.2
            distance = atom_a1 - atom_a2
            resid = atom_a2.get_parent().id[1]

            if resid not in contacts_drug and distance < threshold:
                contacts_drug.append(resid)
                #print(atom_a2.get_id(),atom_a2.get_parent().id[1])
                break
    contacts_drug = sorted(contacts_drug)
    return contacts_drug

def find_chain_mol(pdb, mol_id):
    parser = PDB.PDBParser()
    structure = parser.get_structure("structure", pdb)
    
    for model in structure:
        for chain in model:
            for atom in chain.get_atoms():
                atom_parent = atom.get_parent().resname
                if atom_parent == mol_id:
                    return chain.id
                
            raise ValueError("Molecule not found in the PDB.")

pdb_file = sys.argv[1]
drug_name = sys.argv[2]
drug_chain = find_chain_mol(pdb_file, drug_name)
drug_contacts = contact_res(pdb_file, drug_chain, drug_name)

outstr=''

for res in drug_contacts:
    outstr += '{}{},'.format(res, drug_chain)

outstr = outstr.strip(',')

print(outstr)    
