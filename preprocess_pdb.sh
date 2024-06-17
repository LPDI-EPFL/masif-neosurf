#!/bin/bash

### USER INPUT #################################################################################################

while [[ $# -gt 0 ]]; do
  case $1 in
    -o|--outdir)
      OUTDIR=$2
      shift # past argument
      shift # past value
      ;;
    -l|--ligand)
      LIGAND="$2"
      shift # past argument
      shift # past value
      ;;
    -s|--sdf)
      SDFFILE="$(realpath $2)"
      shift # past argument
      shift # past value
      ;;
    -*|--*)
      echo "Unknown option $1"
      exit 1
      ;;
    *)
      POSITIONAL_ARGS+=("$1") # save positional arg
      shift # past argument
      ;;
  esac
done

set -- "${POSITIONAL_ARGS[@]}" # restore positional parameters
PDBFILE="$(realpath $1)"
NAME_CHAIN=$2


if [ -z $PDBFILE ]; then
    echo "[ERROR] Please provide an input PDB file"
    exit 1
fi

if [ -z $NAME_CHAIN ]; then
    echo "[ERROR] Please provide the protein definition (as <PDBID_CHAINID>)"
    exit 1
fi

if [ -z $OUTDIR ]; then
    echo "[ERROR] Please provide an output directory"
    exit 1
fi

if [ -z $LIGAND ]; then
    echo "[INFO] No small molecule provided"
fi

################################################################################################################

# get directory where this script is located
BASEDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )


# set some variables
MASIF_SOURCE=$BASEDIR/masif/source/
MASIF_SEED_SOURCE=$BASEDIR/masif_seed_search/source/
MASIF_DATA=$BASEDIR/masif/data/
MASIF_TARGETS_DIR=$BASEDIR/masif_seed_search/data/masif_targets

export PYTHONPATH=$PYTHONPATH:$MASIF_SOURCE
export PYTHONPATH=$PYTHONPATH:$OUTDIR
OUTDIR="$(realpath $OUTDIR)"
export TMPDIR=$OUTDIR/tmp/


# Move to output directory
mkdir -p $OUTDIR
cd $OUTDIR

# create directories
mkdir -p data_preparation/00-raw_pdbs/
mkdir -p $TMPDIR

# Link required folders
ln -sf $MASIF_TARGETS_DIR/nn_models
#ln -sf $MASIF_DATA/masif_site/nn_models masif_site_models
#ln -sf $MASIF_DATA/masif_ppi_search/nn_models masif_search_models

# ./data_prepare_one_spec.sh $PDBFILE $NAME_CHAIN $LIGAND $MOL2FILE
echo "Precomputing features on $PDBFILE"
PPI_PAIR_ID=$NAME_CHAIN
PDB_ID=$(echo $NAME_CHAIN| cut -d"_" -f1)
CHAIN1=$(echo $NAME_CHAIN| cut -d"_" -f2)
cp $PDBFILE data_preparation/00-raw_pdbs/$PDB_ID\.pdb
python -W ignore $MASIF_SOURCE/data_preparation/01-pdb_extract_and_triangulate.py $PDB_ID\_$CHAIN1 $LIGAND $SDFFILE
return_code=$?


# Run MaSIF
if [ $return_code -eq 0 ]; then
    python $MASIF_SOURCE/data_preparation/04-masif_precompute.py masif_site $PPI_PAIR_ID
    return_code=$?
fi
if [ $return_code -eq 0 ]; then
    python $MASIF_SOURCE/data_preparation/04-masif_precompute.py masif_ppi_search $PPI_PAIR_ID
    return_code=$?
fi

# ./predict_site.sh $NAME_CHAIN
if [ $return_code -eq 0 ]; then
    echo "Running masif site on $PDBFILE"
    python -W ignore $MASIF_SOURCE/masif_site/masif_site_predict.py nn_models.all_feat_3l.custom_params $NAME_CHAIN
    return_code=$?
fi

# ./color_site.sh $NAME_CHAIN
if [ $return_code -eq 0 ]; then
    export PYTHONPATH=$PYTHONPATH:$MASIF_DATA/masif_site/
    python -W ignore $MASIF_SOURCE/masif_site/masif_site_label_surface.py nn_models.all_feat_3l.custom_params $NAME_CHAIN
    return_code=$?
fi

# ./compute_descriptors.sh $NAME_CHAIN
if [ $return_code -eq 0 ]; then
    echo "Computing descriptors"
    export PYTHONPATH=$PYTHONPATH:$MASIF_DATA/masif_ppi_search/
    python $MASIF_SOURCE/masif_ppi_search/masif_ppi_search_comp_desc.py nn_models.sc05.all_feat.custom_params $NAME_CHAIN
    return_code=$?
fi

# copy files required for seed search
if [ $return_code -eq 0 ]; then
    echo "Creating running directory targets/$NAME_CHAIN "
    mkdir -p targets
    cp -r $MASIF_TARGETS_DIR/targets/template/ targets/$NAME_CHAIN
fi

# return to previous directory
cd -
