import os,sys
sys.path.append('..')
from Bio.PDB import * 
from Bio.PDB.DSSP import DSSP
import numpy as np
from scipy.spatial import cKDTree
from global_vars import dssp_path

"""
This script gives an estimation of the percentage of molecular contacts done by hotspot residue found in beta sheet.
The goal is to discard sheet-based binding seeds with too many contacts done by loop region as they will be cropped
during the downstream seed grafting process.
"""

def contacts_per_dssp(pdb: str, seed_chain: str, target_chain: str, dist_cutoff: float=5.0) -> list:
    """Assign DSSP to seed residues and check how many interactions each residues
    contributes to the overall interface within `dist_cutoff`.

    Args:
        pdb (str): Path to PDB.
        seed_chain (str): Seed chain ID.
        target_chain (str): Target chain ID
        dist_cutoff (float, optional): Cutoff distance in Angstrom within which residues are considered to make contacts. Defaults to 5.0.

    Returns:
        list: List containing annotation for each seed residue in the from ['DSSP', 'ResidueID', 'NumberOfContacts']
    """
    
    parser = PDBParser()
    results = []

    # Load target and seed structures
    try:
        pdb_name = pdb
        pdb_struct = parser.get_structure(pdb_name, pdb_name)
    except:
        print('Error with {}'.format(target_name))

    model_pdb = pdb_struct[0]

    atoms_target = Selection.unfold_entities(model_pdb[target_chain], 'A')
    atoms_seed = Selection.unfold_entities(model_pdb[seed_chain], 'A')
    res_seed = []
    for res in Selection.unfold_entities(model_pdb[seed_chain], 'R'):
        # Only keep canonical amino acids
        try:
            if Polypeptide.three_to_index(res.get_resname()) <= 19:
                res_seed.append(res)
        except KeyError as err:
            print(err, res)
            continue

    dssp_pdb = DSSP(model_pdb, pdb_name, dssp=dssp_path)

    coords_target = [x.get_coord() for x in atoms_target]
    coords_seed = [x.get_coord() for x in atoms_seed]

    ckd = cKDTree(coords_target)
    dists_seed_to_target, r = ckd.query(coords_seed)

    # Get the residues in the interface
    interface = np.where(dists_seed_to_target < dist_cutoff)[0]
    resid_interface = [atoms_seed[x].get_parent().get_id()[1] for x in interface]

    chain_dssp = {}
    chain_dssp_ids = []
    res_seed_ids = [res.get_id() for res in res_seed]
    for key in dssp_pdb.keys():
        if key[0] == seed_chain:# and key[1] in res_A_ids:
            chain_dssp.update({key: dssp_pdb[key]})
            chain_dssp_ids.append(key[1])
    chain_dssp_string = (''.join([i[2] for i in chain_dssp.values()]))
    # Check how many contacts each residue makes and which DSSP type that residue has
    for ix, elem in enumerate(chain_dssp.keys()):
        resid = res_seed[ix].get_id()[1]
        if resid in resid_interface and elem[0] == seed_chain:
            results.append([chain_dssp_string[ix], resid, resid_interface.count(resid)])
    return results

def find_pdb_files():

    list_of_files=[]
    cmd = 'find ./out/* -name "*.pdb" -print '
    for file in os.popen(cmd).readlines():
        name = file[:-1]
        list_of_files.append(name)

    return list_of_files

list_of_files = find_pdb_files()
i=0

with open('dssp.csv','w') as outf:
    outf.write('design,ratio\n')
    for filename in list_of_files:
        dssp_results=contacts_per_dssp(filename, 'B', 'A', 2.5)
        total_contacts=0
        sheet_contacts=0
        for j in dssp_results:
            total_contacts+=j[2]
            if j[0]=='E':
                sheet_contacts+=j[2]
        short=filename.split('/')[-1].strip('.pdb')
        i+=1
        ratio=sheet_contacts/total_contacts
        line=short+','+str(ratio)+'\n'
        outf.write(line)
        print('File '+str(i)+' out of '+str(len(list_of_files))+' done.')
print("Done.")
