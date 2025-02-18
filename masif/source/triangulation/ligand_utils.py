import shutil
from pathlib import Path
from io import StringIO, BytesIO
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.AllChem import AssignBondOrdersFromTemplate
from openbabel import openbabel
import prody
prody.confProDy(verbosity='none')

masif_dir = Path(__file__).parent.parent.parent.resolve()
ligand_expo = np.load(Path(masif_dir, 'data', 'masif_neosurf', 'pdb_ligand_expo.npy'), allow_pickle=True).item()


def neutralize_atoms(mol):
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    # mol_out = Chem.AddHs(mol, addCoords=True)
    mol_out = Chem.Mol(mol)
    return mol_out


def amide_to_single_bond(mol2_outfile):
    mol2_new = []
    bond_record = False
    with open(mol2_outfile, 'r') as f:
        for line in f.readlines():
            line = line.rstrip('\n')

            if line.startswith('@<TRIPOS>'):
                # new data record
                if line.startswith('@<TRIPOS>BOND'):
                    bond_record = True
                else:
                    bond_record = False
                mol2_new.append(line)
                continue

            if bond_record:
                # format: bond_id origin_atom_id target_atom_id bond_type
                bond_id, origin_atom_id, target_atom_id, bond_type = line.split()
                bond_type = '1' if bond_type == 'am' else bond_type
                line = '\t'.join([bond_id, origin_atom_id, target_atom_id, bond_type])

            mol2_new.append(line)

    with open(mol2_outfile, 'w') as f:
        f.write("\n".join(mol2_new))


def extract_ligand(pdb_file, ligand_name, ligand_chain, mol2_outfile, sdf_template=None, patched_mol2_file=None):
    pdb = prody.parsePDB(pdb_file)
    ligand = pdb.select(f'chain {ligand_chain} and resname {ligand_name[:3]}')
    assert ligand is not None and len(ligand) > 0, "Ligand not found"

    out = StringIO()
    prody.writePDBStream(out, ligand)
    rdmol = AllChem.MolFromPDBBlock(out.getvalue(), sanitize=True, removeHs=False)

    try:
        if sdf_template is not None:
            template = Chem.SDMolSupplier(sdf_template)[0]
            template = Chem.AddHs(template)
            rdmol = AllChem.AssignBondOrdersFromTemplate(template, rdmol)
            print(f"[INFO] Inferred ligand connectivity from the provided SDF file")
            
        elif ligand_name in ligand_expo:
            print("Extracting ligand connectivity from PDB Ligand Expo")
            smiles, expo_name = ligand_expo[ligand_name]
            template = AllChem.MolFromSmiles(smiles)
            template = Chem.AddHs(template)
            rdmol = AllChem.AssignBondOrdersFromTemplate(template, rdmol)
            print(f"[INFO] Inferred ligand connectivity from the PDB Ligand Expo (name: {expo_name})")

    except ValueError:
        print("Mismatch between PDB and template ligands. Determining bond types with OpenBabel...")
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("pdb", "sdf")
        obmol = openbabel.OBMol()
        obConversion.ReadString(obmol, out.getvalue())
        sdf_string = obConversion.WriteString(obmol)
        sdf_stream = BytesIO(sdf_string.encode('utf-8'))
        template = list(Chem.ForwardSDMolSupplier(sdf_stream, sanitize=True, removeHs=False))[0]
        rdmol = AllChem.AssignBondOrdersFromTemplate(template, rdmol)
        print(f"[INFO] Inferred ligand connectivity with OpenBabel")

    if patched_mol2_file is not None:
        print("[WARNING] Patched mol2 file provided. It is preferred to use the automatic PDB-to-mol2 conversion. "
              "This option should only be used in cases where the default option fails or yields inconsistent results.")
        shutil.copyfile(patched_mol2_file, mol2_outfile)
    else:
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("pdb", "mol2")
        # obConversion.AddOption("a...")  # read options (preceded by 'a')
        # obConversion.AddOption("xl")  # write options (preceded by 'x')
        ob_mol = openbabel.OBMol()
        obConversion.ReadString(ob_mol, out.getvalue())
        obConversion.WriteFile(ob_mol, mol2_outfile)

    # remove amide bond type because it is not supported by PDB2PQR
    amide_to_single_bond(mol2_outfile)

    rdmol = neutralize_atoms(rdmol)
    assert rdmol.GetNumHeavyAtoms() < rdmol.GetNumAtoms(), \
        print("The molecule must be protonated!")
    
    return rdmol
