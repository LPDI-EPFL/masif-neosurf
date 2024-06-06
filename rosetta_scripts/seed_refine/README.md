# Seed refinement protocol

Although our seed database is large (402M patches), the seeds might not be fully compatible with their respective target. Therefore, it is recommended to proceed to a refinement protocol that combines relax and design with Rosetta.

## Step 0 : Set up the global variable

If not done, set up the global variables (Rosetta path, DSSP path, etc...) and install the necessary python working environment. More information in the [README](../README.md) in the upper folder.

## Step 1 : Preparing the input files

1. Place all your seed PDB in the input folder
2. Run the preparation script by typing the path to the target PDB `-t` and the name of the small molecule `-m`

The script will process the input files and, notably, preparing the small molecule PARAMS file which is necessary to run ligand on RosettaScript.

Example:
```
python3 ./1-prepare.py -t ./1LRY_A.pdb -m BB2
```
**Note:** Some "exotic" atoms have an unpredictable behavior with Rosetta (e.g. CL). Check the formatting of the PARAMS file (e.g. tabulation) and that the atom type is present in the Rosetta radii database (`path_to_rosetta/main/database/scoring/score_functions/sc/sc_radii.lib`)

## Step 2 : Running the seed refinement

Run the seed refinement script with the second command:
```
bash ./2-run.sh
```
Seed refinement typically takes about 10-15 minutes with parallelization depending on the size of your protein complex. 

## Step 3 : Post-processing

Run the post-processing script in order to compute useful metrics like the number of atoms in contact with the small molecule, the similarity score and the percentage of beta sheet contact (for beta sheet-absed seeds only).

Run the command by specifying the seed type and the small molecule name:
```
bash ./3-postprocess.sh S BB2
```

The *similarity score* aims to discard seeds that are highly similar and/or converged to a identical solution. This is based on pairwise sequence alignement and will give a score between 0 (unique) and 1 (identical with another seed).

The *beta sheet contact percentage* is useful to discard seeds whose contacts with the target protein are mostly made by loop regions. Due to their difficulty to be grafted to a scaffold protein, loop regions are cropped in the seed grafting protocol, therefore justifying the need of focusing on seeds that primarily rely on contacts made by beta sheets.

## Step 4 : Analysis and selection

A Jupyter notebook is provided to analyze the Rosetta metrics and post-processing scores. This script will guide you to make a selection of seeds to be grafted. It's recommended to aim for 50-100 helical seeds and 100-200 beta sheet seeds. 
