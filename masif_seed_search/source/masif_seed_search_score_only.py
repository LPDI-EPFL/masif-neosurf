#!/usr/bin/env python

import os
from argparse import ArgumentParser
import shutil
import importlib
from pdb import set_trace

import numpy as np
from scipy.spatial import cKDTree
import pymesh
from Bio.PDB import PDBParser

from alignment_evaluation_nn import AlignmentEvaluationNN
from alignment_utils import get_target_vix, load_protein_pcd, get_patch_coords, compute_nn_score, get_patch_geo
from default_config.masif_opts import masif_opts


if __name__ == "__main__":

    # Inputs
    parser = ArgumentParser()
    parser.add_argument("-p", "--params_file", type=str, required=True)
    parser.add_argument("-t", "--target", type=str, required=True)
    parser.add_argument("-b", "--binder", type=str, required=True)
    args = parser.parse_args()

    custom_params_fn = args.params_file
    target_name = args.target
    binder_name = args.binder

    # Parameters 
    custom_params_obj = importlib.import_module(custom_params_fn, package=None)
    params = custom_params_obj.params
    # from score_params import params

    params['target_surf_dir'] = os.path.join(params['masif_target_root'], masif_opts['ply_chain_dir'])
    params['target_iface_dir'] = os.path.join(params['masif_target_root'], masif_opts['site']['out_pred_dir'])
    params['target_ply_iface_dir'] = os.path.join(params['masif_target_root'], masif_opts['site']['out_surf_dir'])
    params['target_pdb_dir'] = os.path.join(params['masif_target_root'], masif_opts['pdb_chain_dir'])
    params['target_desc_dir'] = os.path.join(params['masif_target_root'], masif_opts['ppi_search']['desc_dir'])
    params['target_precomp_dir'] = os.path.join(params['masif_target_root'], masif_opts['ppi_search']['masif_precomputation_dir'])

    params['binder_surf_dir'] = os.path.join(params['top_seed_dir'], masif_opts['ply_chain_dir'])
    params['binder_iface_dir'] = os.path.join(params['top_seed_dir'], masif_opts['site']['out_pred_dir'])
    params['binder_desc_dir'] = os.path.join(params['top_seed_dir'], masif_opts['ppi_search']['desc_dir'])
    params['binder_precomp_dir'] = os.path.join(params['top_seed_dir'], masif_opts['ppi_search']['masif_precomputation_dir'])
    params['binder_ply_iface_dir'] = os.path.join(params['top_seed_dir'], masif_opts['site']['out_surf_dir'])
    params['binder_pdb_dir'] = os.path.join(params['top_seed_dir'], masif_opts['pdb_chain_dir'])


    # Initialize scoring neural network
    nn_score_atomic = AlignmentEvaluationNN(params['nn_score_atomic_fn'], selected_features=[0,1,2,3], max_npoints=params['max_npoints'])
    nn_score_atomic.restore_model()


    # Load target data
    target_ppi_pair_id = target_name
    target_pid = 'p1'
    target_chain_ix = 1

    target_ply_fn = os.path.join(params['target_ply_iface_dir'], target_name+'.ply')
    mymesh = pymesh.load_mesh(target_ply_fn)
    iface = mymesh.get_attribute('vertex_iface')
    target_coord = get_patch_coords(params['target_precomp_dir'], target_ppi_pair_id, target_pid)

    # Open the pdb structure of the target, load into point clouds for fast access.
    parser = PDBParser()
    target_pdb_path = os.path.join(params['target_pdb_dir'],'{}.pdb'.format(target_name))
    target_struct = parser.get_structure(target_pdb_path, target_pdb_path)
    target_atom_coords = [atom.get_coord() for atom in target_struct.get_atoms() if not atom.get_name().startswith('H') ]
    target_ca_coords = [atom.get_coord() for atom in target_struct.get_atoms() if atom.get_id() == 'CA']
    # Create kdtree search trees (for fast comparision).
    target_ca_pcd_tree = cKDTree(np.array(target_ca_coords))
    target_pcd_tree = cKDTree(np.array(target_atom_coords))


    # Define target patch
    # coord = np.array(target_point_coord)
    coord = np.array(params['target_point']['coord'])
    target_cutoff = params['target_point']['cutoff']
    # find atom indices close to the target.
    dists = np.sqrt(np.sum(np.square(mymesh.vertices - coord), axis=1))
    neigh_indices = np.where(dists < target_cutoff)[0]
    # Get a target vertex for every target site.
    print(coord)
    print(neigh_indices)
    target_vertices = get_target_vix(target_coord, iface, num_sites=params['num_sites'], selected_vertices=neigh_indices)



    # Load the target point cloud, descriptors, interface and mesh.
    target_paths = {}
    target_paths['surf_dir'] = params['target_surf_dir'] 
    target_paths['iface_dir'] = params['target_iface_dir'] 
    target_paths['desc_dir'] = params['target_desc_dir'] 
    # target_pcd, target_desc, target_iface, target_mesh = load_protein_pcd(target_ppi_pair_id, target_chain_ix, target_paths, flipped_features=True, read_mesh=True)
    target_pcd, target_desc, target_iface = load_protein_pcd(target_ppi_pair_id, target_chain_ix, target_paths, flipped_features=True, read_mesh=False)


    # Load binder data
    binder_pdb_path = os.path.join(params['binder_pdb_dir'],'{}.pdb'.format(binder_name))
    binder_ply_fn = os.path.join(params['binder_ply_iface_dir'], binder_name+'.ply')

    binder_ppi_pair_id = binder_name
    binder_pid = 'p1'
    binder_chain_ix = 1

    binder_paths = {}
    binder_paths['surf_dir'] = params['binder_surf_dir'] 
    binder_paths['iface_dir'] = params['binder_iface_dir'] 
    binder_paths['desc_dir'] = params['binder_desc_dir']
    binder_pcd, binder_desc, binder_iface = load_protein_pcd(binder_ppi_pair_id, binder_chain_ix, binder_paths, flipped_features=False, read_mesh=False)


    # Output
    outdir = params['out_dir_template'].format(target_name)
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)

    # Copy the pdb structure and the ply file of the target
    shutil.copy(target_pdb_path, outdir)
    shutil.copy(target_ply_fn, outdir)

    # Copy the pdb structure and the ply file of the binder
    shutil.copy(binder_pdb_path, outdir)
    shutil.copy(binder_ply_fn, outdir)


    # Go through every target site in the target
    for site_ix, site_vix in enumerate(target_vertices):

        score_dict = {'target': target_name, 'binder': binder_name, 'target_vix': site_vix}
        
        site_outdir = os.path.join(outdir, 'site_{}'.format(site_ix))

        if not os.path.exists(site_outdir):
            os.makedirs(site_outdir)
        # Get the geodesic patch and descriptor patch for each target patch
        target_patch, target_patch_descs, target_patch_idx = get_patch_geo(target_pcd, target_coord, site_vix, target_desc, flip_normals=True, outward_shift=params['surface_outward_shift'])

        # Find closest point in binder to define the center of the interface patch on the binder side
        # binder_ply_fn = os.path.join(params['binder_ply_iface_dir'], binder_name+'.ply')
        # binder_mesh = pymesh.load_mesh(binder_ply_fn)
        # dists = np.sqrt(np.sum(np.square(binder_mesh.vertices - mymesh.vertices[site_vix][None, :]), axis=1))
        dists = np.sqrt(np.sum(np.square(np.asarray(binder_pcd.points) - mymesh.vertices[site_vix][None, :]), axis=1))
        binder_center_idx = np.argmin(dists)

        binder_coord =  get_patch_coords(params['binder_precomp_dir'], binder_ppi_pair_id, binder_pid, cv=[binder_center_idx])
        binder_patch, binder_patch_descs, binder_patch_idx = get_patch_geo(binder_pcd, binder_coord, binder_center_idx, binder_desc, outward_shift=params['surface_outward_shift'])

        # Write out the patches themselves
        with open(site_outdir+'/target.xyz', 'w+') as out_patch:
            out_patch.write(f"{len(target_patch.points)}\n\n")
            for point in target_patch.points:
                out_patch.write('H {}, {}, {}\n'.format(point[0], point[1], point[2]))

        with open(site_outdir+'/binder.xyz', 'w+') as out_patch:
            out_patch.write(f"{len(binder_patch.points)}\n\n")
            for point in binder_patch.points:
                out_patch.write('H {}, {}, {}\n'.format(point[0], point[1], point[2]))

        # Compute descriptor distances
        score_dict['center_desc_dist'] = np.sqrt(np.sum(np.square(binder_desc[0][binder_center_idx] - target_desc[0][site_vix])))
        
        # dists = np.sqrt(np.sum(np.square(mymesh.vertices[target_patch_idx][:, None, :] - binder_mesh.vertices[None, :, :]), axis=-1))
        dists = np.sqrt(np.sum(np.square(mymesh.vertices[target_patch_idx][:, None, :] - np.asarray(binder_pcd.points)[None, :, :]), axis=-1))
        correspondence = np.argmin(dists, axis=1) # closest point on binder for each point in target patch
        score_dict['average_patch_desc_dist'] = np.mean(np.sqrt(np.sum(np.square(binder_desc[0][correspondence] - target_desc[0][target_patch_idx]), axis=-1)))

        # Compute neural network score
        target_ckdtree = cKDTree(target_patch.points)  # Make a ckdtree with the target vertices.
        d_vi_at, _ = target_pcd_tree.query(np.asarray(binder_patch.points), k=1)  # Compute the distance between every source_surface_vertices and every target vertex.
        align_scores, point_importance = compute_nn_score(target_patch, binder_patch, None, target_patch_descs, binder_patch_descs, target_ckdtree, nn_score_atomic, d_vi_at, 1.0)
        score_dict['nn_score'] = align_scores[0]
        score_dict['desc_dist_score'] = align_scores[1]


        # Write score file
        score_file = os.path.join(site_outdir, '{}.score'.format(binder_name))
        with open(score_file, 'w+') as out_score:
            for k, v in score_dict.items():
                out_score.write(f"{k}: {v}\n")
