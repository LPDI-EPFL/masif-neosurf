# Running a seed search with _MaSIF-neosurf_

## Preprocess a PDB file

Before we can search for complementary binding sites/seeds, we need to triangulate the molecular surface and compute the initial surface features. The script `preprocess_pdb.sh` takes two required positional arguments: the PDB file and a definition of the chain(s) that will be included.

If a small molecule is part of the molecular surface, we need to tell MaSIF-neosurf where to find it in the PDB file (three letter code + chain) using the `-l` flag. Optionally, we can also provide an SDF file with the `-s` flag that will be used to infer the correct connectivity information (i.e. bond types). This SDF file can be downloaded from the PDB website for example.

Example:
```bash
chmod +x ./preprocess_pdb.sh

# Without connectivitiy information (SDF)

./preprocess_pdb input/6o0k.pdb 6O0K_A -l LBM_A

# With connectivity information (SDF) - Optional

./preprocess_pdb input/6o0k.pdb 6O0K_A -l LBM_A -s input/6O0K_A_LBM.sdf
```

This proccess will ultimately generate the surface descriptors and MaSIF-site interface prediction. The interface prediction can be found in the form of a PLY file in the `./output` folder. More explanations how to read PLY files in PyMOL with our plug-in can be found in `./masif-neosurf/masif/source/`

A working folder for the seed search will be created in the `./targets` directory.

## Performing a seed search

Our work is using a helical and beta sheet-based seed libraries (Provided upon request). Therefore, the parameters and submission files are split in "helical" and "sheet" cases. Running the seed search takes about 10-15 with our seed database.

### Step 1 : Set up your search database

If this is your first run, enter the path to your seed/scaffold/designs database on which the search should be performed. In the parameters file, enter the path to your database in `params['masif_db_root']`. If you don't use our MaSIF-seed database, also change the path within the database in `params['top_seed_dir']`.

### Step 2 : Set up your search parameters

You can either let MaSIF choose the most promising site to search (default) or you can specify a target residue around which a patch fingerprint will be used for the search.
To do so, uncomment the following line and complete the target residue information, namely how far the patch can be from the selected atom (`cutoff`), the residue number in the PDB file (`resid`), the chain it belongs (`chain`) and the name of the target atom in the PDB file (`atom_id`).

```bash
params['target_residue'] = {'cutoff': 3.0, 'resid': 95, 'chain': 'B', 'atom_id': 'OD1'}
```

Secondly, you need to define the cutoffs for the search. Recommended values are provided in the file and can be fine-tuned according to the number of output PDB you aim for.

`params['desc_dist_cutoff']` : Euclidian distance between target and query descriptors
`params['iface_cutoff']` : Interface propensity of the searched patch
`params['nn_score_cutoff']` : Interface post-alignement (IPA) score

### Step 3 : Run the search

The search can be done with CPUs in parallel. In this work, we split the seeds to scan in different lists to speed up the searching process. You can launch the search with the command :

```bash
sbatch run_peptides_[helix/sheet].slurm PDB_ID

#Example

sbatch run_peptides_helix.slurm 6O0K_A
```

You can adapt this script with the job scheduling system used in your HPC environement (Here we use SLURM). For our database of 402M patches, this process can take 5-15 min with parallelization depending on the cutoffs you chose. The matching seed PDBs will be written in an `./out_peptides` folder. The number of output seeds will depend on your cutoffs.

## Post-processing

If you want to reproduce the seed refinement and grafting performed in the publication, please refer to the relevant [README](rosetta_scripts/README.md). For these following steps we usually aim to use 500-1000 helical seeds or 1000-2000 sheet-based seeds. 
