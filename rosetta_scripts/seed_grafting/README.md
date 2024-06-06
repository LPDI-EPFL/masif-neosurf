# Seed grafting protocol

This protocol uses RosettaScript and the MotifGraft mover to graft the seed onto a scaffold protein originating from a database. 

## Step 0 : Set up the environment

If not done, set up the global variables (Rosetta path, scaffold database, DSSP path, etc...) and install the necessary python working environment. More information in the [README](../README.md) in the upper folder.

## Step 1 : Prepare the input files

1. Place your selected seed PDBs in the input folder
2. Run the preparation script by typing the type of seeds `-t` (H for 'Helix' or S for 'Sheet') and the name of the small molecule `-m`

The script will process the input files and, notably, preparing the small molecule PARAMS file which is necessary to run ligand on RosettaScript.

Example:
```
python3 ./1-prepare.py -t H -m BB2
```
**Note:** Some "exotic" atoms have an unpredictable behavior with Rosetta (e.g. CL). Check the formatting of the PARAMS file (e.g. tabulation) and that the atom type is present in the Rosetta radii database (`path_to_rosetta/main/database/scoring/score_functions/sc/sc_radii.lib`)

## Step 2 : Running the seed grafting

Run the seed grafting  script with the second command:
```
bash ./2-run.sh
```
Seed grafting typically takes 30-60 minutes with parallelization depending on the size of your protein target and the number of matching scaffolds for your seed.

## Step 3 : Post-processing

Once the grafting is done, run the post-processing script in order to gather all scores files in one single file with only the necessary metrics. Please, specify the type of seeds `-t` (H for 'Helix' or S for 'Sheet') used for the grafting process:
```
python3 ./3-postprocess.py -t H
```

At this stage, it is also recommended to run a prediction with AlphfaFold2 on the monomeric design in order to select only designs that are properly folded. 

## Step 4 : Analysis and selection

A Jupyter notebook is provided to analyze the Rosetta metrics and to guide you to make a selection of final designs.

## Troubleshooting
If you don't get any grafted designs back, here are a couple possibilities:
1. Your seed might be too long : The longer the seed is the more difficult it is to find a suitable scaffold. Maybe you want to consider some manual cropping to increase the success rate.
2. Your seed topology is not adequate : Some helical seed might be kinked due to a proline for example. Manual cropping can be a solution in these cases.
3. One metric is permanently failing : Check in the output which filter is always failing. Eventually turn off the filtering (`confidence = "0"`) but keep the metrics output to look at them individually.
4. Your small molecule params file has a formatting issue : The generation of the params file can sometimes lead to wrong tabulation with some "exotic" atoms, leading to a failure of the ShapeComplementarity filter for example. Have a look at the warning messages of the sanity check performed during file preparation (`./1-prepare.py`) and solved the tabulation problem manually.
5. The atom names of the small molecule are not matching between the PDBs and params file : Although the preparation script is ensuring compatible names between these two files, check if this problem occurs. Otherwise, bonds will be assigned to the wrong atom coordinates. 
