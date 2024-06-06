import numpy as np
from rdkit import Chem
from rdkit.Chem import Crippen
from rdkit.Chem import BRICS


# Kyte Doolittle scale
kd_scale = {}
kd_scale["ILE"] = 4.5
kd_scale["VAL"] = 4.2
kd_scale["LEU"] = 3.8
kd_scale["PHE"] = 2.8
kd_scale["CYS"] = 2.5
kd_scale["MET"] = 1.9
kd_scale["ALA"] = 1.8
kd_scale["GLY"] = -0.4
kd_scale["THR"] = -0.7
kd_scale["SER"] = -0.8
kd_scale["TRP"] = -0.9
kd_scale["TYR"] = -1.3
kd_scale["PRO"] = -1.6
kd_scale["HIS"] = -3.2
kd_scale["GLU"] = -3.5
kd_scale["GLN"] = -3.5
kd_scale["ASP"] = -3.5
kd_scale["ASN"] = -3.5
kd_scale["LYS"] = -3.9
kd_scale["ARG"] = -4.5


def kd_from_logp(logp, kd_min=-4.5, kd_max=4.5):
    return np.clip(-6.2786 + np.exp(0.4772 * logp + 1.8491), kd_min, kd_max)


def get_fragments(rdmol):
    frags = Chem.GetMolFrags(Chem.FragmentOnBRICSBonds(rdmol), asMols=True, sanitizeFrags=False)

    # find exit atoms to be removed
    def find_exits(rdmol):
        exits = []
        for a in rdmol.GetAtoms():
            if a.GetSymbol() == '*':
                exits.append(a.GetIdx())
        return exits

    out_frags = []
    for frag in frags:
        exits = sorted(find_exits(frag), reverse=True)
        efrag = Chem.EditableMol(frag)
        for idx in exits:
            efrag.RemoveAtom(idx)
        efrag = efrag.GetMol()
        Chem.SanitizeMol(efrag, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
        out_frags.append(efrag)
    return out_frags


def mol_assign_kd(mol):
    atom_name_map = {}
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + 1)  # 0 is reserved as default value
        assert atom.GetAtomMapNum() not in atom_name_map
        # atom_name_map[atom.GetAtomMapNum()] = atom.GetProp('_TriposAtomName')
        atom_name_map[atom.GetAtomMapNum()] = atom.GetPDBResidueInfo().GetName().strip()

    assert 0 not in atom_name_map, print("Index 0 is reserved as default value")

    fragments = get_fragments(mol)
    atom_kd = {}
    for frag in fragments:
        frag_logp = Chem.Crippen.MolLogP(frag)
        frag_kd = kd_from_logp(frag_logp)
        for atom in frag.GetAtoms():
            # atom_name = atom.GetProp('_TriposAtomName')
            atom_mapnum = atom.GetAtomMapNum()
            if atom_mapnum == 0:  # atom doesn't exist in original molecule
                print(atom.GetSymbol())
                continue
            atom_name = atom_name_map[atom_mapnum]
            assert atom_name not in atom_kd
            atom_kd[atom_name] = frag_kd

    return atom_kd


# For each vertex in names, compute
def computeHydrophobicity(names, ligand_code=None, rdmol=None):

    if rdmol is not None:
        mol_kd = mol_assign_kd(rdmol)

    hp = np.zeros(len(names))
    for ix, name in enumerate(names):

        aa = name.split("_")[3]

        if rdmol is not None and aa == ligand_code:
            atom_name = name.split("_")[4]
            hp[ix] = mol_kd[atom_name]
        else:
            hp[ix] = kd_scale[aa]

    return hp

