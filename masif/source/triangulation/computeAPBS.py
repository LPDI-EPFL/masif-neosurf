import os
import numpy
from IPython.core.debugger import set_trace
from subprocess import Popen, PIPE
import pymesh

from default_config.global_vars import apbs_bin, pdb2pqr_bin, multivalue_bin
import random

"""
computeAPBS.py: Wrapper function to compute the Poisson Boltzmann electrostatics for a surface using APBS.
Pablo Gainza - LPDI STI EPFL 2019
This file is part of MaSIF.
Released under an Apache License 2.0
"""

def computeAPBS(vertices, pdb_file, tmp_file_base, mol2_file=None):
    """
        Calls APBS, pdb2pqr, and multivalue and returns the charges per vertex
    """
    fields = tmp_file_base.split("/")[0:-1]
    directory = "/".join(fields) + "/"
    filename_base = tmp_file_base.split("/")[-1]
    pdbname = pdb_file.split("/")[-1]
    args = [
        pdb2pqr_bin, pdbname, filename_base,
        "--ff=PARSE",
        "--whitespace",
        "--noopt",
        "--apbs-input={}.in".format(filename_base),
    ]
    if mol2_file is not None:
        args.append("--ligand=" + mol2_file)
    p2 = Popen(args, stdout=PIPE, stderr=PIPE, cwd=directory)
    stdout, stderr = p2.communicate()

    print("### PDB2PQR ###\n", stderr.decode('utf-8'))
    # from pdb import set_trace; set_trace()

    args = [apbs_bin, filename_base + ".in"]
    p2 = Popen(args, stdout=PIPE, stderr=PIPE, cwd=directory)
    stdout, stderr = p2.communicate()
    vertfile = open(directory + "/" + filename_base + ".csv", "w")
    for vert in vertices:
        vertfile.write("{},{},{}\n".format(vert[0], vert[1], vert[2]))
    vertfile.close()

    print("### APBS ###\n", stderr.decode('utf-8'))
    # from pdb import set_trace; set_trace()

    args = [
        multivalue_bin,
        filename_base + ".csv",
        filename_base + ".dx",
        filename_base + "_out.csv",
    ]
    p2 = Popen(args, stdout=PIPE, stderr=PIPE, cwd=directory)
    stdout, stderr = p2.communicate()

    print("### MULTIVALUE ###\n", stderr.decode('utf-8'))
    # from pdb import set_trace; set_trace()

    # Read the charge file
    chargefile = open(tmp_file_base + "_out.csv")
    charges = numpy.array([0.0] * len(vertices))
    for ix, line in enumerate(chargefile.readlines()):
        charges[ix] = float(line.split(",")[3])

    remove_fn = os.path.join(directory, filename_base)
    os.remove(remove_fn)
    os.remove(remove_fn+'.csv')
    os.remove(remove_fn+'.dx')
    os.remove(remove_fn+'.in')
    os.remove(remove_fn+'_out.csv')

    return charges
