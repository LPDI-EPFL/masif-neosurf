all_pdbs = []
with open ('./seed_benchmark_pdbs.txt') as f: 
    for line in f.readlines(): 
        all_pdbs.append(line.rstrip().upper())

with open('./training.txt') as f: 
    with open('./training_seed_benchmark.txt', 'w+') as of:
        for line in f.readlines(): 
            pdbid = line.split('_')[0]
            if pdbid not in all_pdbs: 
                of.write(line)
            else:
                print(line.rstrip())
