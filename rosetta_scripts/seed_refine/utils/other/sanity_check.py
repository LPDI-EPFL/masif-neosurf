import os, sys

def sanity_check(mol_id):
    error = False
    with open('./prep/{}.params'.format(mol_id)) as inf:
        for line in inf:
            split = line.split()
            if split[0] == 'ATOM':
                if line[4:6] != '  ':
                    print('WARNING: The following tabulation problem has been found in ./prep/{}.params file'.format(mol_id))
                    print(line.strip())
                    print('     ^')
                    error = True
            elif split[0] == 'BOND_TYPE':
                if line[9:11] != '  ':
                    print('WARNING: The following tabulation problem has been found in ./prep/{}.params file'.format(mol_id))
                    print(line.strip())
                    print('          ^')
                    error = True
            elif split[0] == 'ICOOR_INTERNAL':
                if line[14:18] != '    ':
                    print('WARNING: The following tabulation problem has been found in ./prep/{}.params file'.format(mol_id))
                    print(line.strip())
                    print('                 ^')
                    error = True
    if error:
        print('WARNING: Some formatting issues have been identified and can lead to unexpected behaviour when running RosettaScript.')
    else:
        print('Sanity check passed.')
