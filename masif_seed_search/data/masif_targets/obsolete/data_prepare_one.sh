#!/bin/bash
masif_root=../../../
masif_source=$masif_root/masif/source/
masif_seed_source=$masif_root/masif_seed_search/source/
export PYTHONPATH=$masif_source
export masif_matlab
#if [ "$1" == "--file" ]
#then
#   echo "Running masif site on $2"
#   PPI_PAIR_ID=$3
#   PDB_ID=$(echo $PPI_PAIR_ID| cut -d"_" -f1)
#   CHAIN1=$(echo $PPI_PAIR_ID| cut -d"_" -f2)
#   CHAIN2=$(echo $PPI_PAIR_ID| cut -d"_" -f3)
#   FILENAME=$2
#        mkdir -p data_preparation/00-raw_pdbs/
#   cp $FILENAME data_preparation/00-raw_pdbs/$PDB_ID\.pdb
#else
    PPI_PAIR_ID=$1
    PDB_ID=$(echo $PPI_PAIR_ID| cut -d"_" -f1)
    CHAIN1=$(echo $PPI_PAIR_ID| cut -d"_" -f2)
    CHAIN2=$(echo $PPI_PAIR_ID| cut -d"_" -f3)
    python -W ignore $masif_source/data_preparation/00-pdb_download.py $PPI_PAIR_ID
#fi
python -W ignore $masif_source/data_preparation/01-pdb_extract_and_triangulate.py $PDB_ID\_$CHAIN1 $2 $3
python $masif_source/data_preparation/04-masif_precompute.py masif_site $PPI_PAIR_ID
python $masif_source/data_preparation/04-masif_precompute.py masif_ppi_search $PPI_PAIR_ID
# UNCOMMENT THESE LINES IF YOU WANT D-AMINO ACIDS
# Make and run for d protein
#python -W ignore $masif_seed_source/01d-make_d_target.py $PDB_ID\_$CHAIN1
#python $masif_source/data_preparation/04-masif_precompute.py masif_site d$PPI_PAIR_ID
#python $masif_source/data_preparation/04-masif_precompute.py masif_ppi_search d$PPI_PAIR_ID
