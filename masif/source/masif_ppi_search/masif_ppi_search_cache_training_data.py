# Header variables and parameters.
import sys
import pymesh
import os
import numpy as np
from IPython.core.debugger import set_trace
import importlib

from scipy.spatial import cKDTree

from default_config.masif_opts import masif_opts

"""
masif_ppi_search_cache_training_data.py: Function to cache all the training data for MaSIF-search. 
                This function extract all the positive pairs and a random number of negative surfaces.
                In the future, the number of negative surfaces should be increased.
Pablo Gainza - LPDI STI EPFL 2019
Released under an Apache License 2.0
"""

params = masif_opts['ppi_search']

if len(sys.argv) > 0:
    custom_params_file = sys.argv[1]
    custom_params = importlib.import_module(custom_params_file, package=None)
    custom_params = custom_params.custom_params

    for key in custom_params:
        print('Setting {} to {} '.format(key, custom_params[key]))
        params[key] = custom_params[key]

if 'pids' not in params:
    params['pids'] = ['p1', 'p2']

# Read the positive first
parent_in_dir = params['masif_precomputation_dir']

binder_rho_wrt_center = []
binder_theta_wrt_center = []
binder_input_feat = []
binder_mask = []

pos_rho_wrt_center = []
pos_theta_wrt_center = []
pos_input_feat = []
pos_mask = []

neg_rho_wrt_center = []
neg_theta_wrt_center = []
neg_input_feat = []
neg_mask = []

np.random.seed(0)
training_idx = []
val_idx = []
test_idx = []
pos_names = []
neg_names = []

training_list = [x.rstrip() for x in open(params['training_list']).readlines()]
testing_list = [x.rstrip() for x in open(params['testing_list']).readlines()]

idx_count = 0
for count, ppi_pair_id in enumerate(os.listdir(parent_in_dir)):
    print(f"{count} / {len(os.listdir(parent_in_dir))}")
    if ppi_pair_id not in testing_list and ppi_pair_id not in training_list:
        continue
    in_dir = parent_in_dir + ppi_pair_id + '/'
    print(ppi_pair_id)

    # Read binder and pos.
    train_val = np.random.random()
    # Read binder first, which is p1.
    try:
        labels = np.load(in_dir + 'p1' + '_sc_labels.npy')
        # Take the median of the percentile 25 shape complementarity.
        mylabels = labels[0]
        labels = np.median(mylabels, axis=1)

    except Exception as e:
        print('Could not open ' + in_dir + 'p1' + '_sc_labels.npy: ' + str(e))
        continue

    # Read the corresponding ply files.
    fields = ppi_pair_id.split('_')
    ply_fn1 = masif_opts['ply_file_template'].format(fields[0], fields[1])
    ply_fn2 = masif_opts['ply_file_template'].format(fields[0], fields[2])

    # pos_labels: points > max_sc_filt and >  min_sc_filt.
    pos_labels = np.where((labels < params['max_sc_filt']) & (labels > params['min_sc_filt']))[0]
    K = int(params['pos_surf_accept_probability'] * len(pos_labels))
    if K < 1:
        continue
    l = np.arange(len(pos_labels))
    np.random.shuffle(l)
    l = l[:K]
    l = pos_labels[l]

    v1 = pymesh.load_mesh(ply_fn1).vertices[l]
    v2 = pymesh.load_mesh(ply_fn2).vertices

    # For each point in v1, find the closest point in v2.
    kdt = cKDTree(v2)
    d, r = kdt.query(v1)
    # Contact points: those within a cutoff distance.
    contact_points = np.where(d < params['pos_interface_cutoff'])[0]
    try:
        k1 = l[contact_points]
    except:
        set_trace()
    k2 = r[contact_points]

    # For negatives, get points in v2 far from p1.
    try:
        kdt = cKDTree(v1)
        dneg, rneg = kdt.query(v2)
    except:
        set_trace()
    k_neg2 = np.where(dneg > params['pos_interface_cutoff'])[0]

    assert len(k1) == len(k2)
    n_pos = len(k1)

    pid = 'p1'  # Binder is p1
    for ii in k1:
        pos_names.append('{}_{}_{}'.format(ppi_pair_id, pid, ii))

    rho_wrt_center = np.load(in_dir + pid + '_rho_wrt_center.npy')
    theta_wrt_center = np.load(in_dir + pid + '_theta_wrt_center.npy')
    input_feat = np.load(in_dir + pid + '_input_feat.npy')
    mask = np.load(in_dir + pid + '_mask.npy')

    binder_rho_wrt_center.append(rho_wrt_center[k1])
    binder_theta_wrt_center.append(theta_wrt_center[k1])
    binder_input_feat.append(input_feat[k1])
    binder_mask.append(mask[k1])

    # Read pos, which is p2.
    pid = 'p2'

    # Read as positives those points.
    rho_wrt_center = np.load(in_dir + pid + '_rho_wrt_center.npy')
    theta_wrt_center = np.load(in_dir + pid + '_theta_wrt_center.npy')
    input_feat = np.load(in_dir + pid + '_input_feat.npy')
    mask = np.load(in_dir + pid + '_mask.npy')
    pos_rho_wrt_center.append(rho_wrt_center[k2])
    pos_theta_wrt_center.append(theta_wrt_center[k2])
    pos_input_feat.append(input_feat[k2])
    pos_mask.append(mask[k2])

    # Get a set of negatives from  p2.
    np.random.shuffle(k_neg2)
    k_neg2 = k_neg2[:(len(k2))]
    assert (len(k_neg2) == n_pos)
    neg_rho_wrt_center.append(rho_wrt_center[k_neg2])
    neg_theta_wrt_center.append(theta_wrt_center[k_neg2])
    neg_input_feat.append(input_feat[k_neg2])
    neg_mask.append(mask[k_neg2])
    for ii in k_neg2:
        neg_names.append('{}_{}_{}'.format(ppi_pair_id, pid, ii))

    # Training, validation or test?
    if ppi_pair_id in training_list:
        if train_val <= params['range_val_samples']:
            training_idx = np.append(training_idx, np.arange(idx_count, idx_count + n_pos))
        elif train_val > params['range_val_samples']:
            val_idx = np.append(val_idx, np.arange(idx_count, idx_count + n_pos))
    else:
        test_idx = np.append(test_idx, np.arange(idx_count, idx_count + n_pos))

    idx_count += n_pos

if not os.path.exists(params['cache_dir']):
    os.makedirs(params['cache_dir'])

binder_rho_wrt_center = np.concatenate(binder_rho_wrt_center, axis=0)
binder_theta_wrt_center = np.concatenate(binder_theta_wrt_center, axis=0)
binder_input_feat = np.concatenate(binder_input_feat, axis=0)
binder_mask = np.concatenate(binder_mask, axis=0)

pos_rho_wrt_center = np.concatenate(pos_rho_wrt_center, axis=0)
pos_theta_wrt_center = np.concatenate(pos_theta_wrt_center, axis=0)
pos_input_feat = np.concatenate(pos_input_feat, axis=0)
pos_mask = np.concatenate(pos_mask, axis=0)

neg_rho_wrt_center = np.concatenate(neg_rho_wrt_center, axis=0)
neg_theta_wrt_center = np.concatenate(neg_theta_wrt_center, axis=0)
neg_input_feat = np.concatenate(neg_input_feat, axis=0)
neg_mask = np.concatenate(neg_mask, axis=0)

# Remove inputs with NaNs
num_patches_old = len(binder_input_feat)
print(f"Num. patches before removing NaNs: {num_patches_old}")
not_nan = ~np.isnan(binder_input_feat).any(axis=(1, 2))
not_nan = not_nan & (~np.isnan(pos_input_feat).any(axis=(1, 2)))
not_nan = not_nan & (~np.isnan(neg_input_feat).any(axis=(1, 2)))
print(f"Found NaNs in {(~not_nan).sum()} patches.")

pos_names = [x for i, x in enumerate(pos_names) if not_nan[i]]
neg_names = [x for i, x in enumerate(neg_names) if not_nan[i]]
binder_rho_wrt_center = binder_rho_wrt_center[not_nan, :]
binder_theta_wrt_center = binder_theta_wrt_center[not_nan, :]
binder_input_feat = binder_input_feat[not_nan, ...]
binder_mask = binder_mask[not_nan, :]
pos_rho_wrt_center = pos_rho_wrt_center[not_nan, :]
pos_theta_wrt_center = pos_theta_wrt_center[not_nan, :]
pos_input_feat = pos_input_feat[not_nan, ...]
pos_mask = pos_mask[not_nan, :]
neg_rho_wrt_center = neg_rho_wrt_center[not_nan, :]
neg_theta_wrt_center = neg_theta_wrt_center[not_nan, :]
neg_input_feat = neg_input_feat[not_nan, ...]
neg_mask = neg_mask[not_nan, :]
new_idx_to_old_idx = np.arange(num_patches_old)[not_nan]
old_idx_to_new_idx = -1 * np.ones(num_patches_old)
# for new_idx, old_idx in enumerate(new_idx_to_old_idx):
#     old_idx_to_new_idx[old_idx] = new_idx
old_idx_to_new_idx[new_idx_to_old_idx] = np.arange(len(new_idx_to_old_idx))
print("Num. -1 in old_idx_to_new_idx:", (old_idx_to_new_idx == -1).sum())
training_idx = old_idx_to_new_idx[training_idx[not_nan[training_idx.astype(int)]].astype(int)]
val_idx = old_idx_to_new_idx[val_idx[not_nan[val_idx.astype(int)]].astype(int)]
test_idx = old_idx_to_new_idx[test_idx[not_nan[test_idx.astype(int)]].astype(int)]
assert (old_idx_to_new_idx == -1).sum() == (~not_nan).sum()
assert not (training_idx == -1).any()
assert not (val_idx == -1).any()
assert not (test_idx == -1).any()
assert (set(training_idx) | set(val_idx) | set(test_idx)) == set(np.arange(len(binder_mask)))
print(f"Num. patches after removing NaNs: {len(binder_input_feat)}")

np.save(params['cache_dir'] + '/pos_names.npy', pos_names)
np.save(params['cache_dir'] + '/neg_names.npy', neg_names)

print("Read {} negative shapes".format(len(neg_rho_wrt_center)))
print("Read {} positive shapes".format(len(pos_rho_wrt_center)))
np.save(params['cache_dir'] + '/binder_rho_wrt_center.npy', binder_rho_wrt_center)
np.save(params['cache_dir'] + '/binder_theta_wrt_center.npy', binder_theta_wrt_center)
np.save(params['cache_dir'] + '/binder_input_feat.npy', binder_input_feat)
np.save(params['cache_dir'] + '/binder_mask.npy', binder_mask)

np.save(params['cache_dir'] + '/pos_training_idx.npy', training_idx)
np.save(params['cache_dir'] + '/pos_val_idx.npy', val_idx)
np.save(params['cache_dir'] + '/pos_test_idx.npy', test_idx)
np.save(params['cache_dir'] + '/pos_rho_wrt_center.npy', pos_rho_wrt_center)
np.save(params['cache_dir'] + '/pos_theta_wrt_center.npy', pos_theta_wrt_center)
np.save(params['cache_dir'] + '/pos_input_feat.npy', pos_input_feat)
np.save(params['cache_dir'] + '/pos_mask.npy', pos_mask)

np.save(params['cache_dir'] + '/neg_training_idx.npy', training_idx)
np.save(params['cache_dir'] + '/neg_val_idx.npy', val_idx)
np.save(params['cache_dir'] + '/neg_test_idx.npy', test_idx)
np.save(params['cache_dir'] + '/neg_rho_wrt_center.npy', neg_rho_wrt_center)
np.save(params['cache_dir'] + '/neg_theta_wrt_center.npy', neg_theta_wrt_center)
np.save(params['cache_dir'] + '/neg_input_feat.npy', neg_input_feat)
np.save(params['cache_dir'] + '/neg_mask.npy', neg_mask)
# np.save(params['cache_dir']+'/neg_names.npy', neg_names)
