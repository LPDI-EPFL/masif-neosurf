masif_seed_root=$(git rev-parse --show-toplevel)
masif_seed_search_root=$masif_seed_root/masif_seed_search
masif_root=$masif_seed_root/masif
masif_target_root=$masif_seed_search_root/data/masif_targets/
export masif_db_root=../../../../../masif/
masif_source=$masif_root/source/
masif_data=$masif_root/data/
export masif_root
export masif_target_root
export PYTHONPATH=$PYTHONPATH:$masif_source:`pwd`
if [[ $# -eq 0 ]] ; then
    echo 'You must pass as an argument the PDB id of the target and its chain (e.g. 3R2X_AB)'
    exit 0
fi
python -W ignore $masif_seed_search_root/source/masif_seed_search_nn.py params_peptides $1 
