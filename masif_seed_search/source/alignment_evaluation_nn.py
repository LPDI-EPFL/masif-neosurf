import tensorflow as tf
import numpy as np
from pathlib import Path
from scipy.spatial import cKDTree
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score
from tensorflow import keras
import os
import time
#import pandas as pd
import pickle
import sys
from tensorflow.keras.callbacks import Callback

class SaveAtEpoch(Callback):
    def __init__(self, filepath):
        self.filepath = filepath
        self.best_loss = float('inf')
    def on_epoch_end(self, epoch, logs={}):
        self.model.save(self.filepath)
        return



class RocCallback(Callback):
    def __init__(self,training_data,validation_data):
        self.x = training_data[0]
        self.y = training_data[1]
        self.x_val = validation_data[0]
        self.y_val = validation_data[1]

    def on_train_begin(self, logs={}):
        return

    def on_train_end(self, logs={}):
        return

    def on_epoch_begin(self, epoch, logs={}):
        return

    def on_epoch_end(self, epoch, logs={}):
        y_pred_train = self.model.predict_proba(self.x)
        roc_train = roc_auc_score(self.y, y_pred_train[:,1])
        y_pred_val = self.model.predict_proba(self.x_val)
        roc_val = roc_auc_score(self.y_val, y_pred_val[:,1])
        print('\rroc-auc_train: %s - roc-auc_val: %s' % (str(round(roc_train,4)),str(round(roc_val,4))),end=100*' '+'\n')

        # Print the number of mean number of clashes for positives and negatives. 
        positives = np.where(self.y == 1)[0]
        negatives = np.where(self.y == 0)[0]
        print(self.x.shape)
        print("Mean CA clashes positives {}".format(np.mean(self.x[positives,0,-2])))
        print("Mean CA clashes negatives {}".format(np.mean(self.x[negatives,0,-2])))
        print("Mean clashes positives {}".format(np.mean(self.x[positives,0,-1])))
        print("Mean clashes negatives {}".format(np.mean(self.x[negatives,0,-1])))
        return

    def on_batch_begin(self, batch, logs={}):
        return

    def on_batch_end(self, batch, logs={}):
        return

class AlignmentEvaluationNN:
    def __init__(self, model_name, selected_features, max_npoints, n_negatives=100, n_positives=1):
        self.max_npoints = max_npoints
        self.model_name = model_name
        self.n_negatives = n_negatives
        self.n_positives = n_positives
        self.selected_features = selected_features
        config = tf.ConfigProto()
        config.gpu_options.allow_growth = True
        self.session = tf.Session(config=config)

        np.random.seed(42)
        tf.random.set_random_seed(42)

        reg = keras.regularizers.l2(l=0.0)
        model = keras.models.Sequential()
        model.add(keras.layers.Conv1D(filters=16,kernel_size=1,strides=1))
        model.add(keras.layers.BatchNormalization())
        model.add(keras.layers.ReLU())
        model.add(keras.layers.Conv1D(filters=32,kernel_size=1,strides=1))
        model.add(keras.layers.BatchNormalization())
        model.add(keras.layers.ReLU())
        model.add(keras.layers.Conv1D(filters=64,kernel_size=1,strides=1))
        model.add(keras.layers.BatchNormalization())
        model.add(keras.layers.ReLU())
        model.add(keras.layers.Conv1D(filters=128,kernel_size=1,strides=1))
        model.add(keras.layers.BatchNormalization())
        model.add(keras.layers.ReLU())
        model.add(keras.layers.Conv1D(filters=256,kernel_size=1,strides=1))
        model.add(keras.layers.BatchNormalization())
        model.add(keras.layers.ReLU())
        model.add(keras.layers.GlobalAveragePooling1D())
        model.add(keras.layers.Dense(128,activation=tf.nn.relu,kernel_regularizer=reg))
        model.add(keras.layers.Dense(64,activation=tf.nn.relu,kernel_regularizer=reg))
        model.add(keras.layers.Dense(32,activation=tf.nn.relu,kernel_regularizer=reg))
        model.add(keras.layers.Dense(16,activation=tf.nn.relu,kernel_regularizer=reg))
        model.add(keras.layers.Dense(8,activation=tf.nn.relu,kernel_regularizer=reg))
        model.add(keras.layers.Dense(4,activation=tf.nn.relu,kernel_regularizer=reg))
        model.add(keras.layers.Dense(2, activation='softmax'))

        opt = keras.optimizers.Adam(lr=1e-4)
        model.compile(optimizer=opt,loss='sparse_categorical_crossentropy',metrics=['accuracy'])
        self.model = model

    def train_nn(self, train_features, train_labels, test_features, test_labels):

        roc = RocCallback(training_data=(train_features[:,:,self.selected_features], train_labels),
                          validation_data=(test_features[:,:,self.selected_features], test_labels))

#        save = SaveAtEpoch(self.model_name)
        callbacks = [
            keras.callbacks.ModelCheckpoint(filepath=self.model_name,save_best_only=True,monitor='val_loss',save_weights_only=True),
#            keras.callbacks.ModelCheckpoint(filepath=self.model_name,save_best_only=True,monitor='val_loss', save_weights_only=False),
            roc
        ]
        history = self.model.fit(train_features[:,:,self.selected_features],train_labels,batch_size=32,epochs=50,validation_split=0.1,shuffle=True,class_weight={0:1.0/self.n_negatives,1:1.0/self.n_positives}, callbacks=callbacks)

    def restore_model(self):
        self.model.load_weights(self.model_name)
        
    def eval_model(self, features, nn_score_cutoff):
            max_npoints = self.max_npoints

            assert(max_npoints == 100 or max_npoints == 200)
            self.nn_score_cutoff = nn_score_cutoff
            features = np.expand_dims(features, 0)
            n_features = features.shape[2]
            distance = features[0,:,0]
    
            features_trimmed = np.zeros((1,max_npoints,n_features))
            for j,f in enumerate(features):
                if f.shape[0]<=max_npoints:
                    features_trimmed[j,:f.shape[0],:] = f
                else:
                    selected_rows = np.random.choice(f.shape[0],max_npoints,replace=False)
                    features_trimmed[j,:,:] = f[selected_rows]
        
            y_test_pred = self.model.predict(features_trimmed)
            y_test_pred = y_test_pred[:,1].reshape((-1,1))

            point_importance = np.zeros(len(distance))

            # Compute point importance.
            if len(distance) <= max_npoints and y_test_pred[0,0] > self.nn_score_cutoff:
                # Evaluate point by point. 
                for i in range(len(distance)):
                    feat_copy = np.copy(features_trimmed)
                    feat_copy[0,i,:] = 0.0
                    point_val = self.model.predict(feat_copy)
                    point_val = point_val[:,1].reshape((-1,1))
                    point_importance[i] = point_val[0,0] - y_test_pred[0,0]
                # Normalize
                d = point_importance
                const = np.max(np.abs(d))/np.std(d)
                d_std = d/(const*np.std(d))
                point_importance = d_std

            return y_test_pred, point_importance

