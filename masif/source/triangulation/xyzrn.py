from Bio.PDB import *
from default_config.chemistry import radii, polarHydrogens
from IPython.core.debugger import set_trace

"""
xyzrn.py: Read a pdb file and output it is in xyzrn for use in MSMS
Pablo Gainza - LPDI STI EPFL 2019
This file is part of MaSIF.
Released under an Apache License 2.0
"""

def output_pdb_as_xyzrn(pdbfilename, xyzrnfilename, keep_hetatms=None):
    """
        pdbfilename: input pdb filename
        xyzrnfilename: output in xyzrn format.
        keep_hetatms: list of hetatms to keep
    """
    if keep_hetatms is None:
        keep_hetatms = []
    parser = PDBParser()
    struct = parser.get_structure(pdbfilename, pdbfilename)
    outfile = open(xyzrnfilename, "w")
    for atom in struct.get_atoms():
        name = atom.get_name()
        residue = atom.get_parent()
        # Ignore hetatms.
        # if residue.get_id()[0] != " " and residue.get_id()[0][-3:] != 'RC8':
        if residue.get_id()[0] != " " and residue.get_id()[0][-3:] not in keep_hetatms:
            continue
        resname = residue.get_resname()
        reskey = residue.get_id()[1]
        chain = residue.get_parent().get_id()
        atomtype = name[0]

        color = "Green"
        coords = None
        if atomtype in radii and (resname in polarHydrogens or resname in keep_hetatms):
            if atomtype == "O":
                color = "Red"
            if atomtype == "N":
                color = "Blue"
            if atomtype == "H":
                if resname in keep_hetatms:
                    pass
                elif name in polarHydrogens[resname]:
                    color = "Blue"  # Polar hydrogens
            coords = "{:.06f} {:.06f} {:.06f}".format(
                atom.get_coord()[0], atom.get_coord()[1], atom.get_coord()[2]
            )
            insertion = "x"
            if residue.get_id()[2] != " ":
                insertion = residue.get_id()[2]
            full_id = "{}_{:d}_{}_{}_{}_{}".format(
                chain, residue.get_id()[1], insertion, resname, name, color
            )
        if coords is not None:
            outfile.write(coords + " " + radii[atomtype] + " 1 " + full_id + "\n")

