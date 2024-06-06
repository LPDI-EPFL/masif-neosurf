masif_root=../../../
masif_seed_root=../../
masif_source=$masif_root/source/
masif_seed_source=$masif_seed_root/source/
masif_matlab=$masif_root/source/matlab_libs/
masif_data=$masif_root/data/
export PYTHONPATH=$PYTHONPATH:$masif_source:$masif_seed_source
python -W ignore -u $masif_seed_source/seed_search_generate_training_data.py $masif_root/data/masif_ppi_search/ 1000 2000 9 training_data_9A/ $1
