import os,sys,string,gzip
import math
from utils import find_length
from utils import find_min_dist
from utils import find_hotspots
from utils import find_nfrag
sys.path.append('..')
from global_vars import ch_target
from global_vars import ch_seed

################################################################################
# Utility script for cropping sheet-based seeds to minimum                     #
# Copyright A. Marchand (2024)                                                 #
################################################################################

seed_id = sys.argv[1]

seed_file = './tmp/seed_{}.pdb'.format(seed_id)
cropseed_file = './tmp/cropseed_{}.pdb'.format(seed_id)
context_file = './tmp/context_{}.pdb'.format(seed_id)

with open(cropseed_file,'w') as outf:
    for j in range (1, find_nfrag(seed_file)+1,1):
        start, len_seed = find_length(seed_file,ch_seed, j)
        hotspots = []
        for i in range(start,start+len_seed,1):
            dist = find_min_dist(context_file,seed_file,i,ch_target)
            print(i, dist)
            if(dist < 2.5):
                hotspots.append(i)
        if(len(hotspots) >= 1):
            with open(seed_file,'r') as in_pdb:
                for line in in_pdb:
                    if(line[0:4] == "ATOM"):
                        resid = int(line[22:26].strip(' '))
                        if((resid >= hotspots[0]) and (resid <= hotspots[-1])):
                            outf.write(line)
            outf.write('TER\n')
            
