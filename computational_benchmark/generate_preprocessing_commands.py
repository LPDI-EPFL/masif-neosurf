import sys
from pathlib import Path


if __name__ == "__main__":

    input_file = sys.argv[1]
    struct_dir = sys.argv[2]
    out_dir = sys.argv[3]
    output_file = sys.argv[4]

    preprocess_script = Path(__file__).absolute().parent.parent / 'preprocess_pdb.sh'

    params_with_ligand = []
    params_without_ligand = []

    with open(input_file, 'r') as f:

        for line in f.readlines():
            line = line.strip()

            if line.startswith("#"):
                continue

            pdb_id, chain_1, chain_2, drug, sdf_name = line.split(',')

            params_with_ligand.append( (pdb_id.lower(), f"{pdb_id}_{chain_1}", drug, sdf_name) )
            params_with_ligand.append( (pdb_id.lower(), f"{pdb_id}_{chain_2}", drug, sdf_name) )

            params_without_ligand.append( (pdb_id.lower(), f"{pdb_id}_{chain_1}") )
            params_without_ligand.append( (pdb_id.lower(), f"{pdb_id}_{chain_2}") )


    all_commands = ''

    # Process PDBs with ligands
    for p in params_with_ligand:
        _pdb, _chain, _drug, _sdf = p
        preprocess_cmd = f"{preprocess_script} {Path(struct_dir, _pdb).with_suffix('.pdb')} {_chain} -l {_drug} -s {Path(struct_dir, _sdf).with_suffix('.sdf')} -o {Path(out_dir, 'with_ligand')}"
        all_commands += preprocess_cmd + '\n'

    # Process PDBs without ligands
    for p in params_without_ligand:
        _pdb, _chain = p
        preprocess_cmd = f"{preprocess_script} {Path(struct_dir, _pdb).with_suffix('.pdb')} {_chain} -o {Path(out_dir, 'without_ligand')}"
        all_commands += preprocess_cmd + '\n'

    with open(output_file, 'w') as f:
        f.write(all_commands)
    