#!/bin/bash

#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 6G
#SBATCH --time 1:00:00
#SBATCH --array=0-199
#SBATCH --output=./logs/slurm_%A_%a.txt

### DEBUG
#SLURM_ARRAY_TASK_ID=42
####


input_file=$1
output_dir=$2
pdb_dir=$output_dir/pdbs
mkdir -p $pdb_dir

# Get the specific input parameters for the current job array index
line=$((SLURM_ARRAY_TASK_ID / 2 + 1))
echo $line
ppi=($(sed -n "${line}p" "$input_file"))
echo $ppi

pdbid=$(echo $ppi | cut -d'_' -f1)
echo $pdbid

# Separate the input parameters for even and odd job IDs
if (( SLURM_ARRAY_TASK_ID % 2 == 0 )); then
    # Even job IDs process the first chain
    chain=$(echo $ppi | cut -d'_' -f2)
else
    # Odd job IDs process the second chain
    chain=$(echo $ppi | cut -d'_' -f3)
fi

echo $chain

# Download the PDB
# wget -O $pdb_dir/$pdbid.pdb --no-clobber https://files.rcsb.org/download/$pdbid.pdb
wget -O $pdb_dir/$pdbid.pdb https://files.rcsb.org/download/$pdbid.pdb

# Run pre-processing
echo Processing ${pdbid}_${chain}...
../preprocess_pdb.sh $pdb_dir/$pdbid.pdb ${pdbid}_${chain} -o $output_dir
