from alignment_evaluation_nn import AlignmentEvaluationNN
import numpy as np
from default_config.masif_opts import masif_opts
from IPython.core.debugger import set_trace

nn_score_atomic = AlignmentEvaluationNN('models/weights_seed_benchmark_12A_0123', selected_features=[0,1,2,3], max_npoints=200) ## Slightly slower but more accurate.
nn_score_atomic.restore_model()

mytest_list = open('lists/testing.txt').readlines()
np.random.shuffle(mytest_list)
#mytest_list = ['4ZQK_A_B']
true_pred = []
false_pred= []
for myelem in mytest_list:
    myelem = myelem.rstrip()
    try:
        feats = np.load('training_data_12A_seed_benchmark/{}/features.npy'.format(myelem))
    except:
        continue
    feats = np.array(feats)
    labels = np.load('training_data_12A_seed_benchmark/{}/source_patch_rmsds.npy'.format(myelem))
    true_labels = np.where(labels < 2.0)[0]
    false_labels = np.where(labels >= 2.0)[0]
    for ix, feat in enumerate(feats):
        ypred = nn_score_atomic.eval_model(feat, 0.9)
        if labels[ix] < 2.0:
            true_pred.append(np.squeeze(ypred[0]))
        else:
            false_pred.append(np.squeeze(ypred[0]))
    
true_pred = np.asarray(true_pred)
false_pred = np.asarray(false_pred)
import matplotlib.pyplot as plt
plt.hist(true_pred, density=True, bins=100)
plt.hist(false_pred, alpha=0.75, density=True, bins=100)
plt.savefig('test_network.pdf', type='pdf')

import sklearn.metrics
ytrue = np.concatenate([np.ones_like(true_pred), np.zeros_like(false_pred)])
ypred = np.concatenate([true_pred, false_pred])
roc_auc_score = sklearn.metrics.roc_auc_score(ytrue, ypred)
print("ROC AUC SCORE: {}".format(roc_auc_score))
set_trace()
