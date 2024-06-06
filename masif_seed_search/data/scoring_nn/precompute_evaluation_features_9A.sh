masif_root=../../../
masif_source=$masif_root/source/
masif_seed_source=../../source/
masif_data=$masif_root/data/
export PYTHONPATH=$PYTHONPATH:$masif_source:
python -W ignore $masif_seed_source/precompute_evaluation_features.py training_data_9A/ $1
