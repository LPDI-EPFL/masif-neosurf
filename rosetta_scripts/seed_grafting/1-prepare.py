#!/usr/bin/python
import sys,os
import subprocess
import argparse
sys.path.append('..')
from global_vars import rosetta_path
from global_vars import scaff_path
from global_vars import scaff_list
from openbabel import openbabel
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.AllChem import AssignBondOrdersFromTemplate
from utils.utils import sanity_check

# Parsing arguments

parser = argparse.ArgumentParser(description='Preparing the files for running the seed grafting procedure.')
parser.add_argument('-m', '--mol', type=str, required=True, help='Small molecule three-letter code (Example: BB2)')
parser.add_argument('-t', '--type', type=str, required=True, help='Type of seeds being grafted : Helix (H) or sheet (S)')
args = parser.parse_args()

mol_id = args.mol
seed_type = args.type

if seed_type == 'H':
    seed_type = 'helix'
elif seed_type == 'S':
    seed_type = 'sheet'
else:
    raise ValueError("The seed type must be an helix ('H') or a sheet ('S').")

pdb_files = [file for file in os.listdir('./in') if file.endswith('.pdb')]
pdb_files = sorted(pdb_files)

if not pdb_files:
    raise ValueError("There is no file in the input folder. Cannot proceed.")

# Step 1: Prepare the ligand for Rosetta

print('### Step 1: Prepare the ligand file for Rosetta ###')
print('Extracting ligand information from target PDB file')

with open(f"./in/{pdb_files[0]}",'r') as inf:
    pdb_block=''
    success=False
    for line in inf:
        if len(line.split()) > 1:
            if line.split()[0]=='HETATM' and line.split()[3]==mol_id:
                pdb_block+=line
                success=True
pdb_block+='TER'

if not success:
    raise ValueError(f"No ligand with the identifier {mol_id} found in the first input file")

rdmol = AllChem.MolFromPDBBlock(pdb_block, sanitize=True, removeHs=False)
if not rdmol.GetNumHeavyAtoms() < rdmol.GetNumAtoms():
    raise ValueError("The molecule must be protonated!")

print('Conversion to SDF file')
sdf_outfile='./prep/{}.sdf'.format(mol_id)
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

with open(f"./prep/{mol_id}.params", "r") as file:
    lines = file.readlines()
if lines[-1].startswith("PDB_ROTAMERS"):
    lines.pop()
    with open(f"./prep/{mol_id}.params", "w") as file:
        file.writelines(lines)

atom_list=[]
with open(f"./prep/{mol_id}.pdb", "r") as file:
    for line in file:
        split = line.split()
        if len(split) > 4:
            if split[0] == 'HETATM':
                atom_id = line[12:16]
                atom_list.append(atom_id)

print('Performing sanity check on {}.params'.format(mol_id))
sanity_check(mol_id)

print('### Step 1 : DONE ###')

# Step 2 : Preparing input files

print('### Step 2 : Preparing input files ###')
print('Fixing ligand atom names and renumbering input file...')

if not os.path.exists('./prep/in'):
    subprocess.run('mkdir ./prep/in/', shell=True)

for i in range(1,len(pdb_files)+1,1):
    success = False
    j = 0

    if i%10 == 0:
        print('{} files out of {} done.'.format(i, len(pdb_files)))
    
    filename = pdb_files[i-1]
    new_filename = filename.replace('.pdb', '_S{}.pdb'.format(i))
    
    with open('./in/{}'.format(filename), 'r') as inf:
        with open('./prep/in/{}'.format(new_filename),'w') as outf:
            for line in inf:
                split = line.split()
                if len(split) > 4:
                    if split[0] == 'HETATM' and split[3] == mol_id:
                        success = True
                        new_line = line[0:12] + atom_list[j] + line[16:-1] + '\n'
                        j += 1
                    else:
                        new_line = line
                else:
                    new_line = line
                outf.write(new_line)
    
    if not success:
        raise ValueError(f"No ligand with the identifier {mol_id} found in file {filename}")

print('All input files prepared')
print('Creating list of input seeds')

pdb_files = [file for file in os.listdir('./prep/in/') if file.endswith('.pdb')]
pdb_files = sorted(pdb_files)

with open('./prep/inlist.txt', 'w') as output_file:
    for f in pdb_files:
        output_file.write(f'./prep/in/{f}\n')

print('### Step 2: DONE ###')

# Step 3 : Prepare the submission file

print('### Step 3: Prepare the submission file ###')
print("Creating file...")

seed_nr = len(pdb_files)
with open(r"./src/template_batch_{}.txt".format(seed_type, 'r')) as file:
    template_content = file.read()

slurm_txt = template_content.format(seed_nr = seed_nr)

with open('./prep/batch_submitter.slurm','w') as outf:
    outf.write(slurm_txt)

print('### Step 3: Done ###')

# Step 4 : Prepare the grafting file

print('### Step 4: Prepare the grafting file ###')
print("Creating file...")

with open(r"./src/template_grafting_{}.txt".format(seed_type, 'r')) as file:
    template_content = file.read()

slurm_txt = template_content.format(scaff_list = scaff_list, scaff_path = scaff_path, rosetta_path = rosetta_path, mol_id = mol_id )

with open('./prep/grafting_submitter.slurm','w') as outf:
    outf.write(slurm_txt)

print('### Step 4: DONE ###')
