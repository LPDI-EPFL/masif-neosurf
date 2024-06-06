"""
ligand_charges.py
Arne Schneuing, Evgenia Elizarova - LPDI, EPFL, 2023
This file is part of MaSIF-neosurf.
"""

import os
import numpy as np
from scipy.spatial.distance import cdist
from rdkit.RDPaths import RDDataDir
from rdkit.Chem import ChemicalFeatures
from rdkit import Chem
from rdkit.Chem.Features.FeatDirUtilsRD import *

from triangulation.charges_utils import (
    computeAngleDeviation,
    computeAnglePenalty
)


"""
Overwrite the rdkit implementation of GetAcceptor1FeatVects() because of the issue discussed here:
https://github.com/rdkit/rdkit/issues/3433
"""
def GetAcceptor1FeatVects(conf, featAtoms, scale=1.5):
    """
    Get the direction vectors for Acceptor of type 1

    This is a acceptor with one heavy atom neighbor. There are two possibilities we will
    consider here
    1. The bond to the heavy atom is a single bond e.g. CO
     In this case we don't know the exact direction and we just use the inversion of this bond direction
     and mark this direction as a 'cone'
    2. The bond to the heavy atom is a double bond e.g. C=O
     In this case the we have two possible direction except in some special cases e.g. SO2
     where again we will use bond direction

    ARGUMENTS:
    featAtoms - list of atoms that are part of the feature
    scale - length of the direction vector
    """
    assert len(featAtoms) == 1
    aid = featAtoms[0]
    mol = conf.GetOwningMol()
    nbrs = mol.GetAtomWithIdx(aid).GetNeighbors()

    cpt = conf.GetAtomPosition(aid)

    # find the adjacent heavy atom
    heavyAt = -1
    for nbr in nbrs:
        if nbr.GetAtomicNum() != 1:
          heavyAt = nbr
          break

#     singleBnd = mol.GetBondBetweenAtoms(aid, heavyAt.GetIdx()).GetBondType() > Chem.BondType.SINGLE
    singleBnd = mol.GetBondBetweenAtoms(aid,heavyAt.GetIdx()).GetBondType() == Chem.BondType.SINGLE  # MODIFIED HERE

    # special scale - if the heavy atom is a sulfur (we should proabably check phosphorous as well)
    sulfur = heavyAt.GetAtomicNum() == 16

    if singleBnd or sulfur:
        v1 = conf.GetAtomPosition(heavyAt.GetIdx())
        v1 -= cpt
        v1.Normalize()
        v1 *= (-1.0 * scale)
        v1 += cpt
        return ((cpt, v1), ), 'cone'

    # ok in this case we will assume that
    # heavy atom is sp2 hybridized and the direction vectors (two of them)
    # are in the same plane, we will find this plane by looking for one
    # of the neighbors of the heavy atom
    hvNbrs = heavyAt.GetNeighbors()
    hvNbr = -1
    for nbr in hvNbrs:
        if nbr.GetIdx() != aid:
          hvNbr = nbr
          break

    pt1 = conf.GetAtomPosition(hvNbr.GetIdx())
    v1 = conf.GetAtomPosition(heavyAt.GetIdx())
    pt1 -= v1
    v1 -= cpt
    rotAxis = v1.CrossProduct(pt1)
    rotAxis.Normalize()
    bv1 = ArbAxisRotation(120, rotAxis, v1)
    bv1.Normalize()
    bv1 *= scale
    bv1 += cpt
    bv2 = ArbAxisRotation(-120, rotAxis, v1)
    bv2.Normalize()
    bv2 *= scale
    bv2 += cpt
    return (
        (cpt, bv1),
        (cpt, bv2),
    ), 'linear'


def get_acceptor_type(acceptor_atom):
    """
    see: https://www.rdkit.org/docs/source/rdkit.Chem.Features.FeatDirUtilsRD.html

    Acceptor of type 1: acceptor with one heavy atom neighbor.
    Acceptor of type 2: acceptor with two adjacent heavy atoms.
    Acceptor of type 3: acceptor with three adjacent heavy atoms.
    """
    return sum([a.GetSymbol() != 'H' for a in acceptor_atom.GetNeighbors()])


def get_acceptors(mol):
    fdefFile = os.path.join(RDDataDir, 'BaseFeatures.fdef')
    featFact = ChemicalFeatures.BuildFeatureFactory(fdefFile)
    feats = featFact.GetFeaturesForMol(mol)
    return {f.GetAtomIds()[0] for f in feats if
            f.GetType() == "SingleAtomAcceptor"}


def get_donors(mol):
    fdefFile = os.path.join(RDDataDir, 'BaseFeatures.fdef')
    featFact = ChemicalFeatures.BuildFeatureFactory(fdefFile)
    feats = featFact.GetFeaturesForMol(mol)
    return {f.GetAtomIds()[0] for f in feats if
            f.GetType() == "SingleAtomDonor"}


def get_donor_hydrogens(mol):
    donors = get_donors(mol)
    return {atom.GetIdx() for i in donors for atom in
            mol.GetAtomWithIdx(i).GetNeighbors() if atom.GetSymbol() == 'H'}


def mol_from_pdb(pdb_file, ligand_code):
    with open(pdb_file, 'r') as f:
        ligand_lines = [x for x in f.readlines() if 'HETATM' in x and ligand_code in x]
    assert len(set([int(x.split()[5]) for x in ligand_lines])) == 1, \
        print(f"Ligand {ligand_code} not unique in {pdb_file}")
    return Chem.MolFromPDBBlock(''.join(ligand_lines), sanitize=False,
                                removeHs=False)


def prepare_rdmol(mol):
    # name_to_idx = {a.GetProp('_TriposAtomName'): a.GetIdx()
    #                for a in mol.GetAtoms()}
    name_to_idx = {a.GetPDBResidueInfo().GetName().strip(): a.GetIdx()
                   for a in mol.GetAtoms()}
    donor_hydrogens = get_donor_hydrogens(mol)
    acceptors = get_acceptors(mol)
    return mol, name_to_idx, donor_hydrogens, acceptors


def computeChargeHelperMol(rdmol, atom_idx, donor_hydrogens, acceptors, v):

    # Check if it is a polar hydrogen.
    if atom_idx in donor_hydrogens:
        donor = rdmol.GetAtomWithIdx(atom_idx).GetNeighbors()[0]
        a = rdmol.GetConformer().GetAtomPosition(donor.GetIdx())
        b = rdmol.GetConformer().GetAtomPosition(atom_idx)
        # Donor-H is always 180.0 degrees, = pi
        angle_deviation = computeAngleDeviation(a, b, v, np.pi)
        angle_penalty = computeAnglePenalty(angle_deviation)
        return 1.0 * angle_penalty

    # Check if it is an acceptor oxygen or nitrogen
    elif atom_idx in acceptors:

        acceptor_type = get_acceptor_type(rdmol.GetAtomWithIdx(atom_idx))
        vects, vec_type = globals()[f"GetAcceptor{acceptor_type}FeatVects"](
            rdmol.GetConformer(), [atom_idx], scale=1.0)

        # Note: the ideal direction is often ambiguous, e.g. one of two
        #   options if two lone electron pairs are present or anywhere on a
        #   cone.
        #   We compute the smallest possible angle deviation in these cases.
        if vec_type == 'cone':
            # the ideal direction is anywhere on a cone at 60 degrees from
            # this vector
            sp, ep = vects[0]
            angle_deviation = computeAngleDeviation(ep, sp, v, np.pi / 3)
        else:
            # consider all options and select minimum
            angle_deviation = min([computeAngleDeviation(ep, sp, v, 0)
                                   for sp, ep in vects])

        angle_penalty = computeAnglePenalty(angle_deviation)

        # There's no need to compute the plane penalty separately,
        # out-of-plane directions are already penalized
        plane_penalty = 1.0

        return -1.0 * angle_penalty * plane_penalty

    return 0.0
