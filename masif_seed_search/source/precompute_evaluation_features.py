import os
from IPython.core.debugger import set_trace
# import dask
import sys
import numpy as np
from scipy.spatial import cKDTree
import glob
from default_config.masif_opts import masif_opts

"""
precompute_evaluation_features.py: Precompute the features for the actual training from the 'decoy' transformations.
Freyr Sverrisson and Pablo Gainza - LPDI STI EPFL 2019-2020
Released under an Apache License 2.0
"""

# @dask.delayed
def save_nn(d):
    """ Computes nearest neighbours of points on aligned patches on target patch"""

    aligned_source_patches = np.load(
        d + "/aligned_source_patches.npy", encoding="latin1"
    , allow_pickle=True)
    target_patch = np.load(d + "/target_patch.npy", allow_pickle=True)

    num_source_patches = aligned_source_patches.shape[0]
    all_nn_indxs = []
    all_nn_dists = []
    target_tree = cKDTree(target_patch)
    for i in range(num_source_patches):
        nn_dists, nn_inds = target_tree.query(aligned_source_patches[i], k=1)
        all_nn_indxs.append(nn_inds)
        all_nn_dists.append(nn_dists)

    np.save(d + "/all_nn_indxs.npy", all_nn_indxs)
    np.save(d + "/all_nn_dists.npy", all_nn_dists)

    return True


# @dask.delayed
def preprocess_protein(d):
    """ Precomputes features to train evaluate network on"""
    # Features: 
        # 0: distance between matched points, 
        # 1: descriptor distance
        # 2: normal product
        # 3: vertex_atom_distance
    n_features = 4

    aligned_source_patches = np.load(
        d + "/aligned_source_patches.npy", encoding="latin1"
    , allow_pickle=True)
    aligned_source_patches_descs = np.load(
        d + "/aligned_source_patches_descs_0.npy", encoding="latin1"
    , allow_pickle=True)
    aligned_source_patches_normals = np.load(
        d + "/aligned_source_patches_normals.npy", encoding="latin1"
    , allow_pickle=True)
    aligned_source_patches_to_atom_distance = np.load(
        d + "/aligned_source_patches_dists_to_target_atoms.npy", encoding="latin1"
    , allow_pickle=True)

    target_patch_descs = np.load(d + "/target_patch_descs_0.npy", allow_pickle=True)
    target_patch_normals = np.load(d + "/target_patch_normals.npy", allow_pickle=True)

    all_nn_dists = np.load(d + "/all_nn_dists.npy", encoding="latin1", allow_pickle=True)
    all_nn_indxs = np.load(d + "/all_nn_indxs.npy", encoding="latin1", allow_pickle=True)
    protein_npoints = []
    protein_features = []
    for chosen in range(len(aligned_source_patches)):
        npoints = aligned_source_patches[chosen].shape[0]
        nn_dists = all_nn_dists[chosen]
        nn_inds = all_nn_indxs[chosen]
        desc_dists = np.linalg.norm(
            aligned_source_patches_descs[chosen] - target_patch_descs[nn_inds], axis=1
        )
        normal_dp = np.diag(
            np.dot(
                aligned_source_patches_normals[chosen], target_patch_normals[nn_inds].T
            )
        )
        protein_npoints.append(npoints)
        # Features are 1/dist, 1/desc_dist and normal dot product
        features = np.zeros((npoints, n_features))
        nn_dists[nn_dists < 0.5] = 0.5
        features[:npoints, 0] = 1.0 / nn_dists
        features[:npoints, 1] = 1.0 / desc_dists
        features[:npoints, 2] = normal_dp
        # Chemical features for the source patches.
        dist_to_atom = aligned_source_patches_to_atom_distance[chosen]
        dist_to_atom[dist_to_atom < 0.25] = 0.25 
        dist_to_atom = 1.0/dist_to_atom 
        try:
            features[:npoints, 3] = dist_to_atom
        except:
            set_trace()

        protein_features.append(features)
        

    np.save(d + "/features", protein_features)
    return True


data_dir = sys.argv[1]
ppi_pair_id = sys.argv[2]
data_list = glob.glob(data_dir + ppi_pair_id)
results = []
results2 = []
for d in data_list:
    if not os.path.isdir(d):
        continue
    nn_path = d + "/all_nn_dists.npy"
    #if not os.path.exists(nn_path):
    _ = save_nn(d)
        # results.append(save_nn(d))
    features_path = d + "/features.npy"
    #if not os.path.exists(features_path):
    _ = preprocess_protein(d)
        # results2.append(preprocess_protein(d))
# print(len(results), len(results2))
# results = dask.compute(*results)
# results2 = dask.compute(*results2)
