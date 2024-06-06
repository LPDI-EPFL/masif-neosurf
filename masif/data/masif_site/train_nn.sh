masif_root=$(git rev-parse --show-toplevel)
masif_source=$masif_root/masif/source/
masif_data=$masif_root/masif/data/
export PYTHONPATH=$PYTHONPATH:$masif_source:$masif_data/masif_site/
python3 -W ignore $masif_source/masif_site/masif_site_train.py $1
