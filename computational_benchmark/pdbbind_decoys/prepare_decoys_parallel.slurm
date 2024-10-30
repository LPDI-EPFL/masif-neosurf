#!/bin/bash

#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 6G
#SBATCH --time 01:00:00
#SBATCH --array=1-8879
#SBATCH --output=./logs/slurm_%A_%a.txt


input_file=$1
output_dir=$2
pdb_dir=$output_dir/pdbs
mkdir -p $pdb_dir

# Get the specific input parameters for the current job array index
protein=($(sed -n "${SLURM_ARRAY_TASK_ID}p" "$input_file"))
echo $protein
pdbid=$(echo $protein | cut -d'_' -f1)
echo $pdbid
chain=$(echo $protein | cut -d'_' -f2)
echo $chain

# Download the PDB
wget -O $pdb_dir/$pdbid.pdb --no-clobber https://files.rcsb.org/download/$pdbid.pdb

# Run pre-processing
if [ -e ${output_dir}/descriptors/sc05/all_feat/${pdbid}_${chain}/p1_desc_straight.npy ]; then
    echo "${pdbid}_${chain} already processed!"
else
    echo "Processing ${pdbid}_${chain}..."
    ../preprocess_pdb.sh $pdb_dir/$pdbid.pdb ${pdbid}_${chain} -o $output_dir
fi
