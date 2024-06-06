#!/usr/bin/env python
# coding: utf-8
import sys
from geometry.open3d_import import *
#import ipdb
import numpy as np
import os
from sklearn.manifold import TSNE
from Bio.PDB import *
import copy
import scipy.sparse as spio
from default_config.masif_opts import masif_opts
import sys
from IPython.core.debugger import set_trace
from scipy.spatial import cKDTree
import time
import scipy.spatial 
from alignment_utils import * 

print(sys.argv)
if len(sys.argv) != 7:
    print('Usage: {} data_dir K ransac_iter patch_radius output_dir pdb_list_index'.format(sys.argv[0]))
    print('data_dir: Location of data directory.')
    print('K: Number of descriptors to run')
    sys.exit(1)

data_dir = sys.argv[1]
K=int(sys.argv[2])
ransac_iter = int(sys.argv[3])
# set patch radius fixed at 9A
PATCH_RADIUS = float(sys.argv[4])
if PATCH_RADIUS != 12 or PATCH_RADIUS != 9: 
    print("Currently only supporting patches of radius ~9 (about 100 points) and radius 12 (about 200 points)")
out_base = sys.argv[5]
pdb_list_index = int(sys.argv[6])

surf_dir = os.path.join(data_dir,masif_opts['ply_chain_dir'])
  
desc_dir = os.path.join(data_dir,masif_opts['ppi_search']['desc_dir'])

pdb_dir = os.path.join(data_dir, masif_opts['pdb_chain_dir'])
# precomputation dir to obtain the center of the interface (this is only computed for ppi_search)
precomp_dir_12A = os.path.join(data_dir, masif_opts['ppi_search']['masif_precomputation_dir'])
if PATCH_RADIUS == 12:
    precomp_dir = os.path.join(data_dir, masif_opts['ppi_search']['masif_precomputation_dir'])
elif PATCH_RADIUS == 9:
    precomp_dir = os.path.join(data_dir, masif_opts['site']['masif_precomputation_dir'])
else:
    sys.exit(1)


# In[3]:

benchmark_list = 'lists/testing_seed_benchmark.txt'
pdb_list = open(benchmark_list).readlines()
pdb_list = [x.rstrip() for ix, x in enumerate(pdb_list) if ix % 1000 == pdb_list_index ]

# Read all surfaces. 
all_pc = []
all_desc = []

rand_list = np.copy(pdb_list)
np.random.seed(0)

p2_descriptors_straight = []
p2_point_clouds = []
p2_patch_coords = []
p2_names = []


# Read all of p1, the target. p1 will have flipped descriptors.

all_positive_rmsd = []
# Match all descriptors. 
count_found = 0
all_rankings_desc = []
all_time_global = []
for target_ix,target_pdb in enumerate(rand_list):
    outdir = out_base+'/'+target_pdb+'/'
#    if os.path.exists(outdir):
#        continue
    tic = time.time()
    print(target_pdb)
    target_pdb_id = target_pdb.split('_')[0]
    chains = target_pdb.split('_')[1:]
        
    # Load target descriptors for global matching. 
    try:
        target_desc_sc05 = np.load(os.path.join(desc_dir,target_pdb,'p1_desc_flipped.npy'))
    except:
        print('Failed opening descriptors for {}'.format(target_pdb))
        continue
    target_desc = np.array([target_desc_sc05])

    pdb_chain = '_'.join([target_pdb.split('_')[0], target_pdb.split('_')[1]])

    # Load target point cloud
    target_pc = os.path.join(surf_dir,'{}.ply'.format(target_pdb_id+'_'+chains[0]))
    target_pcd = read_point_cloud(target_pc)

    # Read the point with the highest shape compl.
    sc_labels = np.load(os.path.join(precomp_dir_12A,target_pdb,'p1_sc_labels.npy'))
    center_point = np.argmax(np.median(np.nan_to_num(sc_labels[0]),axis=1))

    # Go through each source descriptor, find the top descriptors, store id+pdb
    tic_global = time.time()
    num_negs = 0
    all_desc_dists = []
    all_pdb_id = []
    all_vix = []
    gt_dists = []

    for source_ix, source_pdb in enumerate(rand_list):

        try:
            source_desc = np.load(os.path.join(desc_dir,source_pdb,'p2_desc_straight.npy'))
        except:
            print("Failed opening descritors for source {}".format(source_pdb))
            continue

        
        desc_dists = np.linalg.norm(source_desc - target_desc[0][center_point],axis=1)
        all_desc_dists.append(desc_dists) 
        all_pdb_id.append([source_pdb]*len(desc_dists))
        all_vix.append(np.arange(len(desc_dists)))
        
            
    all_desc_dists = np.concatenate(all_desc_dists, axis =0)
    all_pdb_id = np.concatenate(all_pdb_id, axis = 0)
    all_vix = np.concatenate(all_vix, axis = 0)
    
    ranking = np.argsort(all_desc_dists)

    # Load target geodesic distances.
    target_coord= get_patch_coords(precomp_dir, target_pdb, 'p1',[center_point])

    # Get the geodesic patch and descriptor patch for the target.
    target_patch, target_patch_descs, idxs = \
             get_patch_geo(target_pcd,target_coord,center_point,\
                     target_desc, flip_normals=True)

    ## Load the structures of the target and the source (to get the ground truth).    
    parser = PDBParser()
    target_struct = parser.get_structure('{}_{}'.format(target_pdb_id,chains[0]), os.path.join(pdb_dir,'{}_{}.pdb'.format(target_pdb_id,chains[0])))
    gt_source_struct = parser.get_structure('{}_{}'.format(target_pdb_id,chains[1]), os.path.join(pdb_dir,'{}_{}.pdb'.format(target_pdb_id,chains[1])))
    # Get coordinates of atoms for the ground truth and target. Ignore protons.
    target_atom_coords = [atom.get_coord() for atom in target_struct.get_atoms() if not atom.get_name().startswith('H')]
    target_ca_coords = [atom.get_coord() for atom in target_struct.get_atoms() if atom.get_id() == 'CA']
    target_atom_coord_pcd = PointCloud()
    target_ca_coord_pcd = PointCloud()
    target_atom_coord_pcd.points = Vector3dVector(np.array(target_atom_coords))
    target_ca_coord_pcd.points = Vector3dVector(np.array(target_ca_coords))
    # Create with cKDTree 
    target_atom_pcd_tree = cKDTree(np.array(target_atom_coords))#KDTreeFlann(target_atom_coord_pcd)
    target_ca_pcd_tree = cKDTree(np.array(target_ca_coords))#KDTreeFlann(target_ca_coord_pcd)
    
    found = False
    myrank_desc = float('inf')
    
    chosen_top = ranking[0:K]
    
    pos_rmsd = []
    time_global = time.time() - tic

    # For each source pdb, for each matched patch store: 
    # Random, and alignment transformation. 
    alignment_transformations = []
    # The aligned source patches themselves. 
    aligned_source_patches = []
    aligned_source_patches_normals = []
    aligned_source_patches_descs_0 = []
    aligned_source_patches_dists_to_target_atoms = []

    # The geodesic distance from the center 

    # Chemical descriptors, no sc filter
    # all feature descriptors, no sc filter 
    # The index of the center. 
    center_index_sources = []
    # The RMSDs, inf if not from same complex
    source_patch_rmsds = []
    # The pdb names 
    source_patch_names = []

    
    # Now that the matches have been done, Go thorugh every source pdb again
    for source_ix, source_pdb in enumerate(rand_list):
        tic = time.time()
        viii = chosen_top[np.where(all_pdb_id[chosen_top] == source_pdb)[0]]

        source_vix = all_vix[viii]
        
        if len(source_vix) == 0:
            print("No vertices selected in source {}".format(source_pdb))
            continue

        # Continue with this pdb.    
        pdb_id = source_pdb.split('_')[0]
        chain = source_pdb.split('_')[2]
        print('Reading: {}'.format(source_pdb))
        source_pcd = read_point_cloud(os.path.join(surf_dir,'{}.ply'.format(pdb_id+'_'+chain)))

        source_desc_sc05 = np.load(os.path.join(desc_dir,source_pdb,'p2_desc_straight.npy'))
        try:
            source_coords= get_patch_coords(precomp_dir,source_pdb, 'p2', cv=source_vix)

            source_desc = np.array([source_desc_sc05])

        except:
            print("Failed loading data for {}".format(source_pdb))
            continue
        
        # Randomly rotate and translate.  
        tic = time.time()
        params = {}
        params['ransac_radius'] = 1.5
        params['ransac_iter'] = 2000
        params['surface_outward_shift'] = 0.25
        all_results, all_source_patch, all_source_descs, all_source_idx \
                = multidock(source_pcd, source_coords, source_desc,source_vix, target_patch, \
                        target_patch_descs, params) 
        toc = time.time()
        print('Ransac time {}'.format((toc - tic)))
        num_negs = num_negs

        for j,res in enumerate(all_results):
            if res.fitness < 1e-8:
                continue
            alignment_transformations.append(res.transformation)
            aligned_source_patches.append(np.asarray(all_source_patch[j].points))
            aligned_source_patches_normals.append(np.asarray(all_source_patch[j].normals))
            aligned_source_patches_descs_0.append(np.asarray(all_source_descs[j][0].data).T)
            center_index_sources.append(source_vix[j])
            source_patch_names.append(source_pdb)

        # If this is the source_pdb, get the ground truth. 
        if source_pdb == target_pdb: 
            
            tic_cc = time.time()
            for j,res in enumerate(all_results):
                if res.fitness < 1e-8:
                    print('Low fitness') 
                    continue
                rmsd = test_alignments(res.transformation, gt_source_struct, target_atom_pcd_tree)
                clashing_ca,clashing, vertex_to_atom_distance =  count_clashes(res.transformation, \
                        np.asarray(all_source_patch[j].points), gt_source_struct, \
                        target_ca_pcd_tree, target_atom_pcd_tree)
                aligned_source_patches_dists_to_target_atoms.append(vertex_to_atom_distance)
                if rmsd < 2.0:
                    rank_val = np.where(chosen_top == viii[j])[0][0]
                    pos_rmsd.append(rmsd)
                    print('Found positive, descriptor rank: {} rmsd: {}'.format(rank_val, rmsd))
                    found = True
                    myrank_desc = min(rank_val, myrank_desc)
                    source_patch_rmsds.append(rmsd)
                else:
                    source_patch_rmsds.append(rmsd)
            toc_cc = time.time()
        else:
            tic_cc = time.time()
            tic_load_struct = time.time()
            source_struct = parser.get_structure('{}_{}'.format(pdb_id,chain), os.path.join(pdb_dir,'{}_{}.pdb'.format(pdb_id,chain)))
            for j,res in enumerate(all_results):
                if res.fitness > 0:
                    clashing_ca,clashing, vertex_to_atom_distance =  count_clashes(res.transformation, np.asarray(all_source_patch[j].points), source_struct, target_ca_pcd_tree, target_atom_pcd_tree)
                    aligned_source_patches_dists_to_target_atoms.append(vertex_to_atom_distance)
                    source_patch_rmsds.append(float('inf'))
            toc_cc = time.time()
    if found: 
        count_found += 1
        all_rankings_desc.append(myrank_desc)
        print("descriptor rank: {}".format(myrank_desc))
    else:
        print('N/D')
        
    all_positive_rmsd.append(pos_rmsd)
    all_time_global.append(time_global)
    # Make out directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    # Save data for this source patch. 
    np.save(os.path.join(outdir,'source_patch_names'), source_patch_names)
    np.save(os.path.join(outdir,'alignment_transformations'), alignment_transformations)
    np.save(os.path.join(outdir,'aligned_source_patches'), np.asarray(aligned_source_patches))
    np.save(os.path.join(outdir,'aligned_source_patches_normals'), np.asarray(aligned_source_patches_normals))
    np.save(os.path.join(outdir,'aligned_source_patches_descs_0'), np.asarray(aligned_source_patches_descs_0))
    np.save(os.path.join(outdir,'aligned_source_patches_dists_to_target_atoms'), np.asarray(aligned_source_patches_dists_to_target_atoms))
    np.save(os.path.join(outdir,'center_index_sources'), center_index_sources)
    np.save(os.path.join(outdir,'source_patch_rmsds'), source_patch_rmsds)
    np.save(os.path.join(outdir,'target_patch'), np.asarray(target_patch.points))
    np.save(os.path.join(outdir,'target_patch_normals'), np.asarray(target_patch.normals))
    np.save(os.path.join(outdir,'target_patch_descs_0'), np.asarray(target_patch_descs[0].data).T)
    np.save(os.path.join(outdir,'target_center_point'), center_point)
    np.save(os.path.join(outdir,'source_center_points'), center_index_sources)
            
print('All alignments took {}'.format(np.sum(all_time_global)))


# In[18]:


all_pos = []
all_neg = []
ranks = []
unranked = 0

rmsds = []
num_success =1


