import os,sys,string
import math
from utils import find_length
from utils import find_min_dist
from utils import find_hotspots
sys.path.append('..')
from global_vars import ch_target
from global_vars import ch_seed

################################################################################
# Utility script for cropping helical seeds to minimum                         #
# Copyright A. Marchand (2024)                                                 #
################################################################################

seed_id = sys.argv[1]

hotspots=find_hotspots('./tmp/context_{}.pdb'.format(seed_id),'./tmp/seed_{}.pdb'.format(seed_id), ch_seed, ch_target)

hotspot_list = hotspots.split(':')
hotspot_list = [int(x) for x in hotspot_list]
first_res = hotspot_list[0]
last_res = hotspot_list[-1]
start_resid, len_seed  = find_length('./tmp/seed_{}.pdb'.format(seed_id), ch_seed)

while (last_res - first_res >= 15):
    if (hotspot_list[1] - hotspot_list[0]) > 2:
        hotspot_list = hotspot_list[1:]
        continue
    elif (hotspot_list[-1] - hotspot_list[-2]) > 2:
        hotspot_list = hotspot_list[:-1]
        continue
    else:
        break

with open('./tmp/seed_{}.pdb'.format(seed_id), 'r') as inf:
    with open('./tmp/cropseed_{}.pdb'.format(seed_id),'w') as outf:
        for line in inf:
            elements = line.split()
            if len(elements) > 1:
                if (elements[0] == 'ATOM') and (int(line[22:26].strip(' ')) >= (start_resid + hotspot_list[0] - 1)) and (int(line[22:26].strip(' ')) <= (start_resid + hotspot_list[-1] -1)):
                    outf.write(line)
        outf.write('TER')


