#!/bin/bash
masif_root=../../../
masif_source=$masif_root/masif/source/
export PYTHONPATH=$masif_source:.
python -W ignore $masif_source/masif_site/masif_site_predict.py nn_models.all_feat_3l.custom_params $1 $2
