#!/usr/bin/python
import sys,os
import subprocess
import argparse
sys.path.append('..')
from global_vars import rosetta_path
from openbabel import openbabel
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.AllChem import AssignBondOrdersFromTemplate
from utils.other.sanity_check import *


parser = argparse.ArgumentParser(description='Preparing the fasta file for running ColabFold')
parser.add_argument('-t', '--target', type=str, required=True, help='Path to the protein target PDB file (Example: ../1LRY_A.pdb)')
parser.add_argument('-m', '--mol', type=str, required=True, help='Small molecule three-letter code (Example: BB2)')
args = parser.parse_args()

target_file = args.target
target_name = target_file.split('/')[-1]
mol_id = args.mol

# Step 1 : Prepare the list of seed to refine

print('### Step 1: Created input seed file list ###')

pdb_files = [file for file in os.listdir('./in') if file.endswith('.pdb')]

if not pdb_files:
    raise ValueError("There is no file in the input folder. Cannot proceed.")

print("Creating file...")
with open('./prep/inlist.txt', 'w') as output_file:
    for f in pdb_files:
        output_file.write(f'./prep/{target_name} ./in/{f}\n')

print('### Step 1: DONE ###')

# Step 2 : Prepare the ligand for Rosetta

print('###Step 2: Prepare the ligand and input files for Rosetta ###')
print('Extracting ligand information from target PDB file')

with open(target_file,'r') as inf:
    pdb_block=''
    success=False
    for line in inf:
        if len(line.split())>4:
            if line.split()[0]=='HETATM' and line.split()[3]==mol_id:
                pdb_block+=line
                success=True
pdb_block+='TER'

if not success:
    raise ValueError(f"No ligand with the identifier {mol_id} found in {target_file}")

rdmol = AllChem.MolFromPDBBlock(pdb_block, sanitize=True, removeHs=False)
if not rdmol.GetNumHeavyAtoms() < rdmol.GetNumAtoms():
    raise ValueError("The molecule must be protonated!")

print('Conversion to SDF file')
sdf_outfile='./prep/'+mol_id+'.sdf'
obConversion = openbabel.OBConversion()
obConversion.SetInAndOutFormats("pdb", "sdf")
ob_mol = openbabel.OBMol()
obConversion.ReadString(ob_mol, pdb_block)
obConversion.WriteFile(ob_mol, sdf_outfile)

print('Conversion to PARAM file')
command=f"{rosetta_path}/main/source/scripts/python/public/molfile_to_params.py -n {mol_id} -p {mol_id} --chain A --clobber --conformers-in-one-file ./prep/{mol_id}.sdf"
process = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
output, error = process.communicate()
print(output.decode())

if len(error.decode())>1:
     raise ValueError(error.decode())

subprocess.run(f"mv {mol_id}* ./prep/", shell=True)

print('Performing sanity check on {}.params'.format(mol_id))
sanity_check(mol_id)

with open(f"./prep/{mol_id}.params", "r") as file:
    lines = file.readlines()
if lines[-1].startswith("PDB_ROTAMERS"):
    lines.pop()
    with open(f"./prep/{mol_id}.params", "w") as file:
        file.writelines(lines)

print('Creating new target file with new atom name')
with open(target_file,'r') as old_target, open(f"./prep/{target_name}", 'w') as new_target:
    pasted=False
    for line in old_target:
        if line.split()[0]=='HETATM' and line.split()[3]==mol_id:
            ligand_pdb='./prep/{}.pdb'.format(mol_id)
            with open(ligand_pdb,'r') as ligand_f:
                new_target.write(ligand_f.read())
            break
        else:
            new_target.write(line)

print('### Step 2: DONE ###')

# Step 3 : Prepare the submission file

print('### Step 3: Prepare the submission file ###')
print("Creating file...")

seed_nr = len(pdb_files)
slurm_txt = f'''#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 6000
#SBATCH --time 00:30:00
#SBATCH --array=1-{seed_nr}
#SBATCH --qos=serial
#SBATCH --output=./exelogs/out/splitfile.%A_%a.out
#SBATCH --error=./exelogs/err/splitfile.%A_%a.err

ROSETTA_BIN={rosetta_path}/main/source/bin/rosetta_scripts.linuxiccrelease

LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p ./prep/inlist.txt)
MYID=$SLURM_ARRAY_TASK_ID
SEED_PATH=$(echo $LINE | cut -d ' ' -f 2)
SEED_NAME=$(basename $SEED_PATH | cut -d "." -f1)
TARGET_PATH=$(echo $LINE | cut -d ' ' -f 1)
TARGET_NAME=$(basename $TARGET_PATH | cut -d "." -f1)
IN_FILE=./relaxed/$TARGET_NAME\_$SEED_NAME.pdb
echo $LINE > ./tmp/tmp_list_$MYID.txt | $LINE
FIX_PRE=$(python3 ./utils/other/get_mol_targetres.py $TARGET_PATH {mol_id})

$ROSETTA_BIN -parser:protocol ./src/relax.xml \\
    -list ./tmp/tmp_list_$MYID.txt \\
    -no_nstruct_label \\
    -out:path:all ./relaxed \\
    -overwrite \\
    -renumber_pdb \\
    -parser:script_vars mol={mol_id} fix=$FIX_PRE\\
    -in:file:extra_res_fa ./prep/{mol_id}.params

rm ./tmp/tmp_list_$MYID.txt
FIX_POST=$(python3 ./utils/other/get_mol_targetres.py $IN_FILE {mol_id})

$ROSETTA_BIN -parser:protocol ./src/interface_design.xml \\
    -s $IN_FILE \\
    -no_nstruct_label \\
    -out:path:all ./out \\
    -overwrite \\
    -renumber_pdb \\
    -parser:script_vars mol={mol_id} fix=$FIX_POST\\
    -holes:dalphaball {rosetta_path}/main/source/external/DAlpahBall/DAlphaBall.gcc \\
    -in:file:extra_res_fa ./prep/{mol_id}.params'''

with open('./prep/batch_submitter.slurm','w') as outf:
    outf.write(slurm_txt)

print('### Step 3: DONE ###')
