# Binder recovery benchmark

This directory contains the code required to re-run our computational benchmark 
in which we aim to recover known binding partners from drug-induced protein 
complexes.
To this end, we split 14 ligand-induced protein complexes into two subunits each 
and search against a database of these 28 proteins plus 200 decoys.
Since the ligand can be considered as a part of each subunit, we have 28 
benchmarking cases. 
The search can be highly parallelised and should only take a couple of minutes.


## Summary of targets
- 6QTL : Caffeine recognizing nanobody
- 7DC8 : Switch Ab Fab and hIL6R in complex with ATP
- 4DRI : FKBP51, Rapamycin and the FRB fragment of mTOR
- 6OB5 : Computationally-designed, modular sense/response system (S3-2D)
- 1A7X : FKBP12-FK1012 COMPLEX
- 6ENG : 43K ATPase domain of Escherichia coli gyrase B in complex with an aminocoumarin
- 6H0F : DDB1-CRBN-pomalidomide complex bound to IKZF1(ZF2)
- 6SJ7 : human DDB1-DDA1-DCAF15 E3 ubiquitin ligase bound to RBM39 and Indisulam
- 7TE8 : CA14-CBD-DB21 ternary complex
- 1TCO : CALCINEURIN A FRAGMENT, CALCINEURIN B, FKBP12 AND THE IMMUNOSUPPRESSANT DRUG FK506 (TACROLIMUS) (Use chain C vs AB)
- 1S9D:  ARF1 IN COMPLEX WITH BREFELDIN A AND A SEC7 DOMAIN (Ignore GDP)
- 4MDK : Cdc34-ubiquitin-CC0651 complex
- 3QEL : NMDA receptor subunit GluN1 and GluN2B in complex with ifenprodil
- 6N4N: Crystal structure of the designed protein DNCR2/danoprevir/NS3a complex

[`benchmark_pdbs.txt`](benchmark_pdbs.txt) contains a summary of all benchmark 
complexes, a definition of the two subunits, and the relevant small molecule.
The format is `PDB ID, chain 1, chain 2, drug, SDF name`.
The drug is specified as `<3_letter_code>_<chain>`.

Each PDB can be fetched from RCSB together with the "Instance coordinates" file 
in SDF format (using the online interface or a command like this: 
```html
https://models.rcsb.org/v1/6qtl/ligand?auth_seq_id=201&label_asym_id=J&encoding=sdf&filename=6qtl_J_CFF.sdf
```

## Decoys

We provide two sets of decoys which are taken from known PPIs without small molecules at the 
interface
- [`masif_ppi_search_benchmark_list.txt`](decoys/masif_ppi_search_benchmark_list.txt),
- [`pdbbind_decoy_list.txt`](pdbbind_decoys/pdbbind_decoy_list.txt).


## Preparing the structures

Prior to the search, we need to preprocess all target complexes as well as the 
database proteins. This includes surface triangulation, feature calculation, 
MaSIF-site interface prediction and computation of MaSIF-search descriptors.
Every subunit of the target complexes is processed individually **and** together 
with the drug.
An already processed version of the dataset is provided on [Zenodo](https://zenodo.org/records/11509001).

### Targets

A helper script generates all the required processing commands for the target 
complexes (to process each subunit with and without the small molecule).
```bash
python generate_preprocessing_commands.py benchmark_pdbs.txt <structure_dir> <processed_dir> preprocessing_commands.sh
```
`<structure_dir>` should contain all the PDBs and SDFs of the benchmark complexes.


All commands can be executed like this:
```bash
source preprocessing_commands.sh
```
The outputs will be written to `<processed_dir>`.

***NOTE:*** Two benchmark complexes (`7TE8`/`P0T_C` and `6ENG`/`BHW_B`) unfortunately required manual intervention 
because the automatic protonation and mol2 conversion produced inconsistent results.
We include manually modified mol2 files in the `mol2_files` folder.
These can be passed to the [processing script](../preprocess_pdb.sh) using the the `-m` flag.
We are happy to provide further instructions or the processed files upon request.


### Decoys

We provide a script that parallelises the downloading and preprocessing of the decoy proteins using slurm:
```bash
sbatch pdbbind_decoys/prepare_decoys_parallel.slurm pdbbind_decoys/pdbbind_decoy_list.txt <decoy_dir>
```
All outputs will be saved in `<decoy_dir>`.


## Running the search
To finally run a search, we first create a parameter file for each target 
(with ligand) in `benchmark_pdbs.txt` using:
```bash
python make_param_files.py <structure_dir> <processed_dir> <decoy_dir> <search_output_dir>
```
This script creates the output folder structure including parameter files in `<search_output_dir>`.
The parameters can be changed by modifying the template in `make_param_files.py` and running it again.

The search & docking step can then be launched in parallel using the newly created parameter files, for example:
```bash
sbatch run_search.slurm <search_output_dir>/masif-neosurf-benchmark/with_ligand/search_params_targets
```
This step has to be repeated for all `search_params_*` folders (4 folders should have been created in the previous step).
When all searches are completed, the hits and docked PDB structures can be found in `<search_output_dir>`.
Some helper functions for further analysis are included in [`analysis_utils.py`](analysis_utils.py).
