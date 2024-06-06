#!/bin/bash

#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 6G
#SBATCH --time 3:00:00
#SBATCH --array=0-27
#SBATCH --output=./logs/slurm_%A_%a.txt

### DEBUG
#SLURM_ARRAY_TASK_ID=0
####

BASEDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
BASEDIR=$(realpath $BASEDIR/..)

# clean the environment
unset PYTHONPATH

input_folder=$1

# Get the list of files in the input folder
files=(${input_folder}/params_*)

# Get the specific input file for the current job array index
input_file=${files[$SLURM_ARRAY_TASK_ID]}
input_file=$(basename $input_file)
input_file="${input_file%.*}"  # remove file extension
echo $input_file

pdb_id=$(basename "$input_file" | cut -d'_' -f2)
chain=$(basename "$input_file" | cut -d'_' -f3)
echo ${pdb_id}_$chain

# Param files are loaded as modules and must therefore be available in the python path
export PYTHONPATH=$PYTHONPATH:$input_folder

# Launch MaSIF seed search
masif_root=$BASEDIR/masif
export masif_root

masif_source=$masif_root/source/
export PYTHONPATH=$PYTHONPATH:$masif_source:`pwd`

masif_seed_search_root=$BASEDIR/masif_seed_search

# RUN
python -W ignore $masif_seed_search_root/source/masif_seed_search_nn.py $input_file ${pdb_id}_$chain
