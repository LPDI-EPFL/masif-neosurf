import os
import sys
import shutil
from pathlib import Path
from argparse import ArgumentParser

import pymesh
from scipy.spatial import cKDTree
import numpy as np
from Bio.PDB import PDBParser

# import MaSIF modules
masif_neosurf_dir = Path(__file__).resolve().parent
sys.path.append(str(Path(masif_neosurf_dir, 'masif', 'source').resolve()))
sys.path.append(str(Path(masif_neosurf_dir, 'masif_seed_search', 'source').resolve()))
from masif.source.default_config.masif_opts import masif_opts
from masif_seed_search.source.alignment_evaluation_nn import AlignmentEvaluationNN
from masif_seed_search.source.alignment_utils import get_patch_coords, load_protein_pcd, get_patch_geo, get_target_vix, match_descriptors, align_protein


def masif_search(params):

    # Decide whether to run complementarity or similarity search
    if params['similarity_mode']:
        params['surface_outward_shift'] = 0.0
        params['allowed_CA_clashes'] = float('inf')
        params['allowed_heavy_atom_clashes'] = float('inf')
        params['nn_score_cutoff'] = 0.0  # neural network score does not work here
        flip_target_features = False
        flip_target_normals = False
        print("Running surface similarity search. (experimental feature)")
    else:
        flip_target_features = True
        flip_target_normals = True
        print("Running surface complementarity search.")

    # # Load target patches.
    target_ppi_pair_id = params['target_name']
    target_pid = 'p1'
    target_chain_ix = 1

    seed_ppi_pair_ids = params['seed_ppi_pair_ids']

    # Initialize two neural networks - one that does not account for atomic clashes (initial filter) and one with clashes. 
    nn_score_atomic = AlignmentEvaluationNN(params['nn_score_atomic_fn'], selected_features=[0,1,2,3], max_npoints=params['max_npoints'])
    nn_score_atomic.restore_model()

    # Go through every 12A patch in the target protein -- get a sorted least in order of the highest iface mean in the patch
    target_ply_fn = os.path.join(params['target_ply_iface_dir'], target_ppi_pair_id + '.ply')
    mymesh = pymesh.load_mesh(target_ply_fn)
    iface = mymesh.get_attribute('vertex_iface')
    target_coord = get_patch_coords(params['target_precomp_dir'], target_ppi_pair_id, target_pid)


    # Define target and source paths (for interface scores, descriptors, ply files)
    target_paths = {}
    target_paths['surf_dir'] = params['target_surf_dir'] 
    target_paths['iface_dir'] = params['target_iface_dir'] 
    target_paths['desc_dir'] = params['target_desc_dir'] 

    source_paths = {}
    source_paths['surf_dir'] = params['seed_surf_dir'] 
    source_paths['iface_dir'] = params['seed_iface_dir'] 
    source_paths['desc_dir'] = params['seed_desc_dir'] 

    # Load the target point cloud, descriptors, interface and mesh.
    target_pcd, target_desc, target_iface, target_mesh = load_protein_pcd(target_ppi_pair_id, target_chain_ix, target_paths, flipped_features=flip_target_features, read_mesh=True)

    # Open the pdb structure of the target, load into point clouds for fast access.
    parser = PDBParser()
    target_pdb_path = os.path.join(params['target_pdb_dir'], '{}.pdb'.format(target_ppi_pair_id))
    target_struct = parser.get_structure(target_pdb_path, target_pdb_path)
    target_atom_coords = [atom.get_coord() for atom in target_struct.get_atoms() if not atom.get_name().startswith('H') ]
    target_ca_coords = [atom.get_coord() for atom in target_struct.get_atoms() if atom.get_id() == 'CA']
    # Create kdtree search trees (for fast comparision).
    target_ca_pcd_tree = cKDTree(np.array(target_ca_coords))
    target_pcd_tree = cKDTree(np.array(target_atom_coords))

    # If a specific residue is selected, then go after that residue
    if 'target_residue' in params:
        target_chain = params['target_residue']['chain']
        target_cutoff = params['target_residue']['cutoff']
        target_atom_id = params['target_residue']['atom_id']
        # Use the tuple for biopython: (' ', resid, ' ')
        target_resid = [x.id for x in target_struct[0][target_chain].get_residues() if x.id[1] == params['target_residue']['resid']]
        assert len(target_resid) == 1, print(f"Target residue ID not unique: {target_resid}")
        target_resid = target_resid[0]
        print(f"Using residue: {target_resid}")
        coord = target_struct[0][target_chain][target_resid][target_atom_id].get_coord()
        # find atom indices close to the target.
        dists = np.sqrt(np.sum(np.square(mymesh.vertices - coord), axis=1))
        neigh_indices = np.where(dists<target_cutoff)[0]
        # Get a target vertex for every target site.
        target_vertices = get_target_vix(target_coord, iface,num_sites=params['num_sites'],selected_vertices=neigh_indices)

    elif 'target_point' in params:
        coord = np.array(params['target_point']['coord'])
        target_cutoff = params['target_point']['cutoff']

        # find atom indices close to the target.
        dists = np.sqrt(np.sum(np.square(mymesh.vertices - coord), axis=1))
        neigh_indices = np.where(dists < target_cutoff)[0]
        # Get a target vertex for every target site.
        target_vertices = get_target_vix(target_coord, iface, num_sites=params['num_sites'], selected_vertices=neigh_indices)

    else:
        # Get a target vertex for every target site.
        target_vertices = get_target_vix(target_coord, iface,num_sites=params['num_sites'])

    outdir = os.path.join(params['out_dir'], params['target_name'])
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)

    # Copy the pdb structure and the ply file of the target
    shutil.copy(target_pdb_path, outdir)
    shutil.copy(target_ply_fn, outdir)

    # Go through every target site in the target
    if 'selected_site_ixs' in params:
        site_ixs = params['selected_site_ixs']
        site_vixs = [vix for ix, vix in enumerate(target_vertices) if ix in site_ixs]
    else:
        site_ixs = [ix for ix, vix in enumerate(target_vertices)]
        site_vixs = target_vertices

    # Go through every selected site
    for site_ix, site_vix in zip(site_ixs,site_vixs):
        site_outdir = os.path.join(outdir, 'site_{}'.format(site_ix))
        if not os.path.exists(site_outdir):
            os.makedirs(site_outdir)

        # Get the geodesic patch and descriptor patch for each target patch
        target_patch, target_patch_descs, target_patch_idx = get_patch_geo(
            target_pcd, target_coord, site_vix, target_desc, 
            flip_normals=flip_target_normals, 
            outward_shift=params['surface_outward_shift'])

        # Make a ckdtree with the target vertices.
        target_ckdtree = cKDTree(target_patch.points)

        # Write out the patch itself.
        with open(site_outdir + '/target.vert', 'w+') as out_patch:
            for point in target_patch.points:
                out_patch.write('{}, {}, {}\n'.format(point[0], point[1], point[2]))

        # Match the top descriptors in the database based on descriptor distance.
        print('Starting to match target descriptor to descriptors from {} proteins; this may take a while.'.format(len(seed_ppi_pair_ids)))
        matched_dict = match_descriptors(seed_ppi_pair_ids, ['p1', 'p2'], target_desc[0][site_vix], params)

        if len(matched_dict.keys())==0:
            continue

        print(" ")
        print("Second stage of MaSIF seed search: each matched descriptor is aligned and scored; this may take a while..")
        count_matched_fragments = 0
        for ix, name in enumerate(matched_dict.keys()):
            align_protein(
                name, \
                target_patch, \
                target_patch_descs, \
                target_ckdtree, \
                target_ca_pcd_tree, \
                target_pcd_tree, \
                source_paths, \
                matched_dict,\
                nn_score_atomic, \
                site_outdir, \
                params
            )
            if (ix + 1) % 1000 == 0:
                print('So far, MaSIF has aligned {} fragments from {} proteins.'.format(count_matched_fragments, ix + 1))
            count_matched_fragments += len(matched_dict[name])


if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument("--database", dest="top_seed_dir", type=Path)
    parser.add_argument("--subset", dest="database_subset", type=Path, default=None)
    parser.add_argument("--out_dir", type=Path)
    parser.add_argument("--target_dir", dest="masif_target_root", type=Path)
    parser.add_argument("--target", dest="target_name", type=str)
    parser.add_argument("--cutoff", dest="target_cutoff", type=float, default=10.0)
    parser.add_argument("--num_sites", type=int, default=1)
    parser.add_argument("--resid", dest="target_resid", type=int, default=None)
    parser.add_argument("--chain", dest="target_chain", type=str, default=None)
    parser.add_argument("--atom_id", dest="target_atom_id", type=str, default=None)
    parser.add_argument("--coord", dest="target_coord", type=float, nargs="+", default=None)
    parser.add_argument("--desc_dist_cutoff", type=float, default=2.0, help="Recommended values: [1.5-2.0] (lower is stricter)")
    parser.add_argument("--iface_cutoff", type=float, default=0.75, help="Recommended values: [0.75-0.95] range (higher is stricter)")
    parser.add_argument("--nn_score_cutoff", type=float, default=0.8, help="# Recommended values: [0.8-0.95] (higher is stricter)")
    parser.add_argument("--desc_dist_score_cutoff", type=float, default=0.0, help="# Recommended values: [0.0-20.0] (higher is stricter)")
    parser.add_argument("--allowed_CA_clashes", type=int, default=0)
    parser.add_argument("--allowed_heavy_atom_clashes", type=int, default=5)
    parser.add_argument("--sim", dest="similarity_mode", action="store_true")
    args = parser.parse_args()

    # Definition of the target patch
    assert ((args.target_resid is not None) and (args.target_chain is not None) and (args.target_atom_id is not None)) ^ (args.target_coord is not None)
    if args.target_coord is not None:
        args.target_point = {'cutoff': args.target_cutoff, 'coord': args.target_coord}
    else:
        args.target_residue = {'cutoff': args.target_cutoff, 'resid': args.target_resid, 'chain': args.target_chain, 'atom_id': args.target_atom_id}


    # Database locations
    args.seed_surf_dir = os.path.join(args.top_seed_dir, masif_opts['ply_chain_dir'])
    args.seed_iface_dir = os.path.join(args.top_seed_dir, masif_opts['site']['out_pred_dir'])
    args.seed_ply_iface_dir = os.path.join(args.top_seed_dir, masif_opts['site']['out_surf_dir'])
    args.seed_pdb_dir = os.path.join(args.top_seed_dir, masif_opts['pdb_chain_dir'])
    args.seed_desc_dir = os.path.join(args.top_seed_dir, masif_opts['ppi_search']['desc_dir'])

    # 9 A
    # args.seed_precomp_dir = os.path.join(args.top_seed_dir,masif_opts['site']['masif_precomputation_dir'])
    # 12 A
    args.seed_precomp_dir = os.path.join(args.top_seed_dir, masif_opts['ppi_search']['masif_precomputation_dir'])

    # Database subset
    if args.database_subset is None:
        # use all available protein IDs in the database
        args.seed_ppi_pair_ids = np.array(os.listdir(args.seed_desc_dir))
    else:
        # use a predefined subset of protein IDs
        with open(args.database_subset, "r") as f:
            args.seed_ppi_pair_ids = [x.rstrip() for x in f.readlines()]

    # Target locations
    args.top_target_dir = os.path.join(args.masif_target_root)
    args.target_surf_dir = os.path.join(args.top_target_dir, masif_opts['ply_chain_dir'])
    args.target_iface_dir = os.path.join(args.masif_target_root, masif_opts['site']['out_pred_dir'])
    args.target_ply_iface_dir = os.path.join(args.masif_target_root, masif_opts['site']['out_surf_dir'])
    args.target_pdb_dir = os.path.join(args.top_target_dir, masif_opts['pdb_chain_dir'])
    args.target_desc_dir = os.path.join(args.top_target_dir, masif_opts['ppi_search']['desc_dir'])
    args.target_desc_dir = os.path.join(args.top_target_dir, masif_opts['ppi_search']['desc_dir'])
    
    # 9 A
    # args.target_precomp_dir = os.path.join(args.top_target_dir,masif_opts['site']['masif_precomputation_dir'])
    # 12 A
    args.target_precomp_dir = os.path.join(args.top_target_dir, masif_opts['ppi_search']['masif_precomputation_dir'])

    # Some hard-coded parameters
    # Neural network scores.
    args.nn_score_atomic_fn = os.path.join(masif_neosurf_dir, "masif_seed_search/data/scoring_nn/models_std/weights_12A_0129.hdf5")
    args.max_npoints = 200

    # Ransac parameters
    args.ransac_iter = 2000
    # Ransac radius - should not be changed.
    args.ransac_radius = 1.5
    # How much to expand the surface for alignment.
    args.surface_outward_shift = 0.25

    masif_search(vars(args))
