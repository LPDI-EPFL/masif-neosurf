import os
from default_config.masif_opts import masif_opts
# CONFIGURATION FILE FOR MASIF_SEED_SEARCH RUN. 

params = {}
# Directory where the database is located.
params['masif_db_root'] = "/work/upcorreia/databases/masif_db_beta_seeds/masif/"
# Seeds (i.e., the fragments) that will be used for this search. 
params['top_seed_dir'] = os.path.join(params['masif_db_root'], 'data/masif_beta_seeds/')
# Root of the targets directory (where to find the sources.)
params['masif_target_root'] = os.environ['masif_target_root']

# Output directory (target_name, target_site, matched_seed)
params['out_dir_template'] = 'out_peptides/{}/'

# Target a specific residue.
# cutoff: the maximum distance of the center of the patch to this residue.
# resid: the residue number in the PDB chain for this residue.
#params['target_residue'] = {'cutoff': 3.0, 'resid': 95, 'chain': 'B', 'atom_id': 'OD1'}


###
# Score cutoffs -- these are empirical values, if they are too loose, then you get a lot of results. 
# Descriptor distance cutoff for the patch. All scores below this value are accepted for further processing.
params['desc_dist_cutoff'] = 2.0 # Recommended values: [1.5-2.0] (lower is stricter)
# Interface cutoff value, all patches with a score below this value are discarded.
params['iface_cutoff'] = 0.8 # Recommended values: [0.75-0.95] range (higher is stricter)
# Neural network score cutoff - Discard anything below this score
params['nn_score_cutoff'] = 0.9 # Recommended values: [0.8-0.95] (higher is stricter)
# Number of clashes to tolerate.
params['allowed_CA_clashes'] = 0
params['allowed_heavy_atom_clashes'] = 5


# Number of sites to target in the protein
params['num_sites'] = 1

# Target a site close to a specific residue

################################
# More advanced configuration (should not be changed.)

# Neural network scores.
params['nn_score_atomic_fn'] = 'models/weights_12A_0129.hdf5'
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
