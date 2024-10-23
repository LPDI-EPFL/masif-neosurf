import sys
from pathlib import Path
import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
import pymesh


basedir = Path(__file__).absolute().parent.parent


global_options = {
    'run_name': 'masif-neosurf-benchmark',
    'desc_dist_cutoff': 3.0,
    'iface_cutoff': 0.0,
    'nn_score_cutoff': 0.8,
    'selection_cutoff': 10.0,
    'num_sites': 3,
    'nn_score_weights': Path(basedir, 'masif_seed_search/data/scoring_nn/models_std/weights_12A_0129.hdf5')
}


template = """
import os
from default_config.masif_opts import masif_opts
# CONFIGURATION FILE FOR MASIF_SEED_SEARCH RUN '{run_name}'. 

params = {{}}
# Directory where the database is located. 
params['top_seed_dir'] = "{database}"
# Root of the targets directory (where to find the sources.)
params['masif_target_root'] = "{target_dir}"

# Output directory (target_name, target_site, matched_seed)
params['out_dir_template'] = "{result_dir}/{{}}/"

# Target a specific residue.
# cutoff: the maximum distance of the center of the patch to this residue.
# resid: the residue number in the PDB chain for this residue.
# params['target_residue'] = {{'cutoff': 10.0, 'resid': {resid}, 'chain': '{chain}', 'atom_id': '{atom_id}'}}
params['target_point'] = {{'cutoff': {selection_cutoff}, 'coord': {coord}}}


###
# Score cutoffs -- these are empirical values, if they are too loose, then you get a lot of results. 
# Descriptor distance cutoff for the patch. All scores below this value are accepted for further processing.
params['desc_dist_cutoff'] = {desc_dist_cutoff} # Recommended values: [1.5-2.0] (lower is stricter)
# Interface cutoff value, all patches with a score below this value are discarded.
params['iface_cutoff'] = {iface_cutoff} # Recommended values: [0.75-0.95] range (higher is stricter)
# Neural network score cutoff - Discard anything below this score
params['nn_score_cutoff'] = {nn_score_cutoff} # Recommended values: [0.8-0.95] (higher is stricter)
# Number of clashes to tolerate.
params['allowed_CA_clashes'] = 0
params['allowed_heavy_atom_clashes'] = 5


# Number of sites to target in the protein
params['num_sites'] = {num_sites}


################################
# More advanced configuration (should not be changed.)

# Neural network scores.
params['nn_score_atomic_fn'] = "{nn_score_weights}"
params['max_npoints'] = 200

# Seed locations
params['seed_surf_dir'] = os.path.join(params['top_seed_dir'], masif_opts['ply_chain_dir'])
params['seed_iface_dir'] = os.path.join(params['top_seed_dir'],masif_opts['site']['out_pred_dir'])
params['seed_ply_iface_dir'] = os.path.join(params['top_seed_dir'],masif_opts['site']['out_surf_dir'])
params['seed_pdb_dir'] = os.path.join(params['top_seed_dir'],masif_opts['pdb_chain_dir'])
params['seed_desc_dir'] = os.path.join(params['top_seed_dir'],masif_opts['ppi_search']['desc_dir'])
# Here is where you set up the radius - right now at 9A.
#params['seed_precomp_dir'] = os.path.join(params['top_seed_dir'],masif_opts['site']['masif_precomputation_dir'])
# 12 A
params['seed_precomp_dir'] = os.path.join(params['top_seed_dir'],masif_opts['ppi_search']['masif_precomputation_dir'])

# Target locations
params['top_target_dir'] = os.path.join(params['masif_target_root'])
params['target_surf_dir'] = os.path.join(params['top_target_dir'], masif_opts['ply_chain_dir'])
params['target_iface_dir'] = os.path.join(params['masif_target_root'],masif_opts['site']['out_pred_dir'])
params['target_ply_iface_dir'] = os.path.join(params['masif_target_root'],masif_opts['site']['out_surf_dir'])
params['target_pdb_dir'] = os.path.join(params['top_target_dir'],masif_opts['pdb_chain_dir'])
params['target_desc_dir'] = os.path.join(params['top_target_dir'],masif_opts['ppi_search']['desc_dir'])
params['target_desc_dir'] = os.path.join(params['top_target_dir'],masif_opts['ppi_search']['desc_dir'])
# 9 A
#params['target_precomp_dir'] = os.path.join(params['top_target_dir'],masif_opts['site']['masif_precomputation_dir'])
# 12 A
params['target_precomp_dir'] = os.path.join(params['top_target_dir'],masif_opts['ppi_search']['masif_precomputation_dir'])

# Ransac parameters
params['ransac_iter'] = 2000
# Ransac radius - should not be changed.
params['ransac_radius'] = 1.5
# How much to expand the surface for alignment.
params['surface_outward_shift'] = 0.25

#params['seed_pdb_list'] = ['3R2X000_C', '3R2X001_C', '3R2X002_C', '3R2X003_C']
"""


if __name__ == "__main__":

    THRESH = 4.0
    patch_def = 'iface_com'  # 'drug_center', 'iface_center_atom', 'iface_com'

    input_file = 'benchmark_pdbs.txt'
    original_structures = sys.argv[1]
    targets_with_ligands = Path(sys.argv[2], 'with_ligand')
    targets_without_ligands = Path(sys.argv[2], 'without_ligand')
    decoys = sys.argv[3]
    output_dir = sys.argv[4]

    output_dir = Path(output_dir, global_options['run_name'])
    out_with_ligand = Path(output_dir, 'with_ligand')
    out_without_ligand = Path(output_dir, 'without_ligand')

    search_setups = [
        {'output_dir': Path(out_with_ligand, 'search_params_targets'), 'database': targets_without_ligands, 'target_dir': targets_with_ligands, 'result_dir': Path(out_with_ligand, 'results_targets')},
        {'output_dir': Path(out_with_ligand, 'search_params_decoys'), 'database': decoys, 'target_dir': targets_with_ligands, 'result_dir': Path(out_with_ligand, 'results_decoys')},
        {'output_dir': Path(out_without_ligand, 'search_params_targets'), 'database': targets_without_ligands, 'target_dir': targets_without_ligands, 'result_dir': Path(out_without_ligand, 'results_targets')},
        {'output_dir': Path(out_without_ligand, 'search_params_decoys'), 'database': decoys, 'target_dir': targets_without_ligands, 'result_dir': Path(out_without_ligand, 'results_decoys')},
    ]

    for i, setup in enumerate(search_setups):
        print(f"\nSetting up search[{i}] with parameters: {setup}")

        with open(input_file, 'r') as f:

            for line in f.readlines():
                line = line.strip()

                if line.startswith("#"):
                    continue

                pdb_id, chain_1, chain_2, drug, sdf_name = line.split(',')

                drug_name, target_chain = drug.split('_')

                pdbfile = Path(original_structures, f"{pdb_id.lower()}.pdb")
                pdb_model = PDBParser(QUIET=True).get_structure('', pdbfile)[0]

                ligand = [res for res in pdb_model[target_chain].get_residues() if res.get_resname() == drug_name]
                assert len(ligand) == 1
                ligand = ligand[0]

                target_resid = ligand.id[1]

                for target_chains, binder_chains in [chain_1, chain_2], [chain_2, chain_1]:
                    print(pdb_id, target_chains, drug)

                    # Get atom furthest away from the protein
                    if patch_def == 'drug_center':
                        protein_atoms = np.stack([a.get_coord() for c in target_chains for a in pdb_model[c].get_atoms()], axis=0)
                        ligand_atoms = np.stack([a.get_coord() for a in ligand.get_atoms()], axis=0)
                        dist = np.sqrt(np.sum((ligand_atoms[:, None, :] - protein_atoms[None, :, :])**2, axis=-1))
                        target_atom_id = np.argmax(np.min(dist, axis=1))
                        target_atom = list(ligand.get_atoms())[target_atom_id]
                        target_atom_name = target_atom.get_name()
                        target_atom_coord = target_atom.get_coord().tolist()

                    # Get atom in the middle of the interface
                    elif patch_def == 'iface_center_atom':
                        target_atoms = [a for c in target_chains for a in pdb_model[c].get_atoms()] + [a for a in ligand.get_atoms()]
                        binder_atoms = [a for c in binder_chains for a in pdb_model[c].get_atoms()]
                        coords_target = np.stack([a.get_coord() for a in target_atoms], axis=0)
                        coords_binder = np.stack([a.get_coord() for a in binder_atoms], axis=0)
                        dist = np.sqrt(np.sum((coords_target[:, None, :] - coords_binder[None, :, :])**2, axis=-1))
                        iface = np.any(dist < THRESH, axis=1)
                        dist_within_iface = np.sqrt(np.sum((coords_target[iface][:, None, :] - coords_target[iface][None, :, :])**2, axis=-1))
                        iface_id = np.argmin(np.max(dist_within_iface, axis=1))
                        target_atom_id = np.arange(len(iface))[iface][iface_id]
                        target_atom = target_atoms[target_atom_id]
                        target_atom_name = target_atom.get_name()
                        target_atom_coord = target_atom.get_coord().tolist()
                        target_resid = target_atom.parent.id[1]
                        target_chain = target_atom.parent.parent.id 

                    # Get atom in the middle of the interface
                    # This is defined based on the original complex which means the point will be identical regardless of the presence of the drug
                    elif patch_def == 'iface_com':
                        # target_atoms = [a for c in target_chains for a in pdb_model[c].get_atoms()] + [a for a in ligand.get_atoms()]
                        # binder_atoms = [a for c in binder_chains for a in pdb_model[c].get_atoms()]
                        target_atoms = [a for c in target_chains for a in pdb_model[c].get_atoms() if is_aa(a.parent)] + [a for a in ligand.get_atoms()]
                        binder_atoms = [a for c in binder_chains for a in pdb_model[c].get_atoms() if is_aa(a.parent)]
                        assert all([is_aa(a.parent) or a.parent.get_resname() == drug_name for a in target_atoms])
                        assert all([is_aa(a.parent) for a in binder_atoms])
                        coords_target = np.stack([a.get_coord() for a in target_atoms], axis=0)
                        coords_binder = np.stack([a.get_coord() for a in binder_atoms], axis=0)
                        dist = np.sqrt(np.sum((coords_target[:, None, :] - coords_binder[None, :, :])**2, axis=-1))
                        iface = np.any(dist < THRESH, axis=1)

                        target_atom_name = 'Coord'
                        target_atom_coord = np.mean(coords_target[iface], axis=0).tolist()
                        target_resid = 'None'
                        target_chain = 'None'
                    
                    else:
                        raise NotImplementedError()

                    # Make sure there are enough surface points in the selections area
                    ply_file = Path(setup['target_dir'], f'output/all_feat_3l/pred_surfaces/{pdb_id}_{target_chains}.ply')
                    surf_mesh = pymesh.load_mesh(str(ply_file))
                    _dist = np.sqrt(np.sum((surf_mesh.vertices - np.array(target_atom_coord)[None, :])**2, axis=-1))
                    within_cutoff = _dist < global_options['selection_cutoff']
                    if np.sum(within_cutoff) < global_options['num_sites']:
                        target_atom_coord = surf_mesh.vertices[np.argmin(_dist)].tolist()
                        print("Not enough surface vertices in selection. Picking closest vertex as center point instead.")

                    options = {
                        'database': setup['database'],
                        'target_dir': setup['target_dir'],
                        'result_dir': setup['result_dir'],
                        'resid': target_resid,
                        'chain': target_chain,
                        'atom_id': target_atom_name,
                        'coord': target_atom_coord,
                        **global_options
                    }

                    Path(setup['output_dir']).mkdir(exist_ok=True, parents=True)
                    output_file = Path(setup['output_dir'], f"params_{pdb_id}_{target_chains}.py")
                    with open(output_file, 'w') as fo:
                        fo.write(template.format(**options)) 

                    # write xyz file of target atom (defining the input patch) for visualization purposes
                    xyz_dir = Path(setup['output_dir'], 'target_points')
                    xyz_dir.mkdir(exist_ok=True)
                    xyz_file = Path(xyz_dir, f"target_point_{pdb_id}_{target_chains}.xyz")
                    with open(xyz_file, 'w') as f:
                        f.write("1\n\n")  # number of 'atoms'
                        f.write(f"X {target_atom_coord[0]:.3f} {target_atom_coord[1]:.3f} {target_atom_coord[2]:.3f}\n")
