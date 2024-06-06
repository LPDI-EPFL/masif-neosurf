masif_root=$(git rev-parse --show-toplevel)
masif_source=$masif_root/masif/source/
masif_data=$masif_root/masif/data/
export PYTHONPATH=$PYTHONPATH:$masif_source:$masif_data/masif_ppi_search/
python3 $masif_source/masif_ppi_search/masif_ppi_search_cache_training_data.py $1
