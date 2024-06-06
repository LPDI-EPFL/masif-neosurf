import os,sys,string
import math
from utils import find_hotspots
sys.path.append('..')
from global_vars import ch_target
from global_vars import ch_seed

################################################################################
# Utility script to output a list of hotspot residues from the seed            #
# Note : The hotspot residues are numbered from 1, not the first residue ID    #
# Copyright A. Marchand (2024)                                                 #
################################################################################

seed_id = sys.argv[1]

hotspots=find_hotspots('./tmp/context_{}.pdb'.format(seed_id),'./tmp/cropseed_{}.pdb'.format(seed_id), ch_seed, ch_target)
print(hotspots)
