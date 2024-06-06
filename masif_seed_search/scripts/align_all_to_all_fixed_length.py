from Bio.PDB import *
from IPython.core.debugger import set_trace
from scipy.spatial import cKDTree
import numpy as np 
import os
import sys

###### METHODS #########
# Save a fragment 
def save_fragment(out_fn, residues):
    structBuild = StructureBuilder.StructureBuilder()
    structBuild.init_structure("output")
    structBuild.init_seg(" ")
    structBuild.init_model(0)
    outputStruct = structBuild.get_structure()
    structBuild.init_chain('A')
    for res in residues: 
        outputStruct[0]['A'].add(res)

    # Output the selected residues
    pdbio = PDBIO()
    pdbio.set_structure(outputStruct)
    pdbio.save(out_fn)

# Return the fragment of residues touching the target dots.
def find_seed_helix(target_dots, residues, frag_size=12):
    ckd = cKDTree(target_dots)
    if len(residues) < frag_size:
        return None
    res_values = []
    for ix, res in enumerate(residues):
        vals_for_res = 0
        for atom in res: 
            dist, r = ckd.query(atom.get_coord())
            if dist < 3.0:
                vals_for_res += 1.0
        res_values.append(np.sum(vals_for_res)) 
    res_values = np.array(res_values)
    best_frag = -1
    best_frag_val = -1
    for ix in range(len(residues) - frag_size +1 ):
        if np.sum(res_values[ix:ix+frag_size]) > best_frag_val:
            best_frag_val = np.sum(res_values[ix:ix+frag_size])
            best_frag = ix
    return residues[best_frag:best_frag+frag_size]


def compute_helix_to_peptide_rmsd(ref_peptide, helix):
    # Compute nearest neighbors between the two peptides and the calpha to calpha rmsd. 

    ref_coord = np.array([x['CA'].get_coord() for x in ref_peptide])
    helix_coord = np.array([x['CA'].get_coord() for x in helix])

    dists = np.sqrt(np.sum(np.square(ref_coord - helix_coord), axis=1))
    rmsd = np.sqrt(np.mean(np.square(dists)))
    return rmsd
###### END METHODS #########

system = sys.argv[1]
orig_pep_index = int(sys.argv[2]) - 1


all_dots = []
with open('out_peptides/{}/site_0/target.vert'.format(system)) as f:
    dots_lines = f.readlines()
    for point_line in dots_lines: 
        point = [float(x) for x in point_line.split(',')]
        all_dots.append(point)

all_dots = np.array(all_dots)

parent_dir = 'out_peptides/{}/site_0'.format(system)

with open('analysis_fixed_size/list_out_peptides.txt') as f:
    list_peptides = f.readlines()
    lp = [x.rstrip() for x in list_peptides]
    list_peptides = lp

current_pep = list_peptides[orig_pep_index]
current_pep_dir = '_'.join(current_pep.split('_')[0:2])

current_struct_fn = os.path.join(parent_dir, current_pep_dir, current_pep+'.pdb')

# Open original peptide.
parser = PDBParser()
orig_struct = parser.get_structure(current_struct_fn, current_struct_fn)
orig_residues = Selection.unfold_entities(orig_struct, 'R')

# Find the amino acids touching target dots 
orig_fragment = find_seed_helix(all_dots, orig_residues)

save_fragment('analysis_fixed_size/out_pdb/{}_frag.pdb'.format(current_pep), orig_fragment)

all_rmsd = []
for mypep in list_peptides:
    mypep_dir = '_'.join(mypep.split('_')[0:2])
    mystruct = os.path.join(parent_dir, mypep_dir, mypep+'.pdb')

    scorefn = mystruct.replace('.pdb', '.score')
    compared_struct = parser.get_structure(mystruct, mystruct)
    residues_compared_struct = Selection.unfold_entities(compared_struct, 'R')

    # Find the amino acids touching target dots 
    compared_fragment = find_seed_helix(all_dots, residues_compared_struct)
    if compared_fragment is None: 
        all_rmsd.append(-1)
        continue
    
    # Do rmsds in both directions and take the minimum.
    rmsd1 = compute_helix_to_peptide_rmsd(orig_fragment, compared_fragment)
    rmsd2 = compute_helix_to_peptide_rmsd(compared_fragment, orig_fragment)

    rmsd = min([rmsd1, rmsd2])
    if rmsd > 100:
        print('ERROR: {} {}'.format(rmsd, num_compared))
        

    all_rmsd.append(rmsd)

np.save('analysis_fixed_size/out_data/rmsd_{}'.format(orig_pep_index), all_rmsd)
