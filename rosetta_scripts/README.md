# Rosetta Scripts
## Requirements
+ Rosetta : Rosetta modelling suite v3.13 onwards. More information [here](https://new.rosettacommons.org/docs/latest/getting_started/Getting-Started#local-installation-and-use-of-rosetta)
+ DSSP : `sudo apt-get install dssp`
+ Scaffold database : See "Data and material availability" statement at the end of the publication to access the different databases.

## Install environment
First, create a working environment by using the YML file provided in this folder:
```
conda env create -f neosurf.yml
conda activate neosurf
```

## Set up global variables
Set up the global variables in the `global_vars.py` file. Please, specify the following elements:
+ Path to Rosetta
+ Path to DSSP binary file
+ Path to the split lists of your scaffold database (See below)
+ The number of split lists that constitute your scaffold database
+ (Optional) Target and seed chain name (*Note:* Chains are renumbered that way during the seed refinement step)

## Creating a scaffold database
Once you downloaded a scaffold database, create some split files that will distribute the path of these scaffold PDBs in different list. This will allow a parallelization of the scaffold search.

1. Create a list with absolute path of all your scaffold PDBs.
2. Run the following command to create some split files:
   ```
   split -l 15 your_final_list.list --additional-suffix .split
   ls -v *.split| cat -n | while read n f; do mv -n "$f" "splitfile_$n"; done
   ```
3. Indicate the path to these split files and adapt the total number of split files in the global variables accordingly.

## Run a seed refinement
For a seed refinement, please read the relevant [README](./seed_refine/README.md)

## Run a seed grafting
After seed refinement and selection of seed candidates, please read the relevant [README](./seed_grafting/README.md) for grafting them onto a scaffold protein. 
