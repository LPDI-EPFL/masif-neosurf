import sys
from utils import find_nfrag

################################################################################
# Utility script to extract output the number of fragments in the seed PDBs    #
# Copyright A. Marchand (2024)                                                 #
################################################################################

seed_id = sys.argv[1]
nfrag = find_nfrag('./tmp/cropseed_{}.pdb'.format(seed_id))

print(nfrag)
