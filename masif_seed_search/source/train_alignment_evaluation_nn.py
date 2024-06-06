import tensorflow as tf
import numpy as np
from pathlib import Path
from scipy.spatial import cKDTree
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score
from tensorflow import keras
from IPython.core.debugger import set_trace
import os
import time
#import pandas as pd
import pickle
import sys
from alignment_evaluation_nn import AlignmentEvaluationNN

n_positives = 1
n_negatives = 100
max_rmsd = 2.0
n_features = 4

def compile_features(data_list):
    all_features = np.empty((len(data_list)*(n_positives+n_negatives),max_npoints,n_features))
    all_labels = np.empty((len(data_list)*(n_positives+n_negatives),1))
    all_npoints = []
    all_idxs = []
    all_nsources = []
    n_samples = 0
    for i,d in enumerate(data_list):
        if (i%100==0) and (i==0):
            print(i,'Feature array size (MB)',all_features.nbytes*1e-6)
            start = time.time()
        elif i%100==0:
            end = time.time()
            print(i,'Feature array size (MB)',all_features.nbytes*1e-6, 'Time', end - start)
            start = time.time()

        
        source_patch_rmsds = np.load(d/'source_patch_rmsds.npy')


        positive_alignments = np.where(source_patch_rmsds<max_rmsd)[0]
        negative_alignments = np.where(source_patch_rmsds>=max_rmsd)[0]

        if len(positive_alignments)==0:#<n_positives:
            continue
        if len(negative_alignments)< n_negatives:#<n_positives:
            continue

        chosen_positives = np.random.choice(positive_alignments,n_positives,replace=False)
        chosen_negatives = np.random.choice(negative_alignments,n_negatives,replace=False)
        chosen_alignments = np.concatenate([chosen_positives,chosen_negatives])
        try:    
            features = np.load(d/'features.npy',encoding='latin1',allow_pickle=True)
        except:
            print("Problem with {}".format(d))
            continue
        n_sources = len(features)
        features = features[chosen_alignments]
        features_trimmed = np.zeros((len(chosen_alignments),max_npoints,n_features))
        for j,f in enumerate(features):
            if f.shape[0]<=max_npoints:
                features_trimmed[j,:f.shape[0],:f.shape[1]] = f
            else:
                selected_rows = np.random.choice(f.shape[0],max_npoints,replace=False)
                features_trimmed[j,:,:f.shape[1]] = f[selected_rows]
            
        labels = np.array((source_patch_rmsds[chosen_alignments]<max_rmsd).astype(int)).reshape(-1,1)
        
        all_features[n_samples:n_samples+len(chosen_alignments),:,:] = features_trimmed
        all_labels[n_samples:n_samples+len(chosen_alignments)] = labels
        n_samples += len(chosen_alignments)


    all_features = all_features[:n_samples]
    all_labels = all_labels[:n_samples]

    all_idxs = np.concatenate([(n_positives+n_negatives)*[i] for i in range(int(all_features.shape[0]/(n_positives+n_negatives)))])
    
    return all_features, all_labels 



selected_features = sys.argv[2].split('_')
selected_features = [int(s) for s in selected_features]


with open('lists/training_seed_benchmark.txt') as f:
    training_list = f.read().splitlines()
    
with open('lists/testing_seed_benchmark.txt') as f:
    testing_list = f.read().splitlines()

patch_size = sys.argv[1]
if patch_size == '12A':
    max_npoints = 200
elif patch_size == '9A':
    max_npoints = 100
else: 
    sys.exit(1)

data_dir = 'training_data_{}_seed_benchmark'.format(patch_size)

data_dir = Path(data_dir)
data_list = list(data_dir.glob('*'))
train_data_list = [d for d in data_list if (os.path.exists(d/'features.npy')) and str(d).split('/')[-1] in training_list]
test_data_list = [d for d in data_list if (os.path.exists(d/'features.npy')) and str(d).split('/')[-1] in testing_list]

train_all_features, train_all_labels = compile_features(train_data_list)
test_all_features, test_all_labels = compile_features(test_data_list)

selected_feat_char = ''.join([str(x) for x in selected_features])
model_name = 'models/weights_seed_benchmark_{}_{}'.format(patch_size,selected_feat_char)

nn = AlignmentEvaluationNN(model_name,selected_features, max_npoints, n_negatives, n_positives)
nn.train_nn(train_all_features, train_all_labels, test_all_features, test_all_labels)


