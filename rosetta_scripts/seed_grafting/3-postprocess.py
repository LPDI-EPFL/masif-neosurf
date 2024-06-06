import os, sys
import pandas as pd
import argparse
import re

# Parsing arguments

parser = argparse.ArgumentParser(description='Gathering all the score files in a unique one.')
parser.add_argument('-t', '--type', type=str, required=True, help='Type of seeds being grafted : Helix (H) or sheet (S)')
args = parser.parse_args()

seed_type = args.type

# Gathering all score files

cmd1 = 'ls ./out | wc -l'
nseed= int(os.popen(cmd1).readlines()[0])
cmd2 = 'ls ./out/seed_1 | wc -l'
narrays = int(os.popen(cmd2).readlines()[0])

df=pd.DataFrame([])

if seed_type == 'H':
    for seedid in range(1,nseed+1,1):
        for array in range(1,narrays+1):
            path='./out/seed_'+str(seedid)+'/'+str(array)+'/score_S'+str(seedid)+'.sc'
            if (os.path.exists(path)):
                df1 = pd.read_csv(path, sep='\s+', skiprows=1).dropna()
                df1 = df1[df1.total_score != 'total_score'].drop(columns=['SCORE:'])
                df1['seed']=str(seedid)
                df1['array']=str(array)
                df = pd.concat([df,df1])
elif seed_type == 'S':
    for seedid in range(1,nseed+1,1):
        for array in range(1,narrays+1):
            path='./out/seed_'+str(seedid)+'/'+str(array)+'/score_S'+str(seedid)+'.sc'
            if (os.path.exists(path)):
                f = open(path,'r')
                outname='tmp_'+str(seedid)+'_'+str(array)+'.sc'
                tmpf = open(outname,'w')
                for line in f:
                    metrics=re.split(' +',line.strip())
                    new_metrics=metrics.copy()
                    if metrics[0]=='SCORE:':
                        for element in metrics:
                            if (element=='0') or (',' in element):
                                new_metrics.remove(element)
                            if (element in ['graft_in_motif_ranges','graft_in_scaffold_ranges', 'graft_out_scaffold_ranges', 'graft_scaffold_size_change']):
                                new_metrics.remove(element)
                    for i in range(0,len(new_metrics),1):
                        tmpf.write(new_metrics[i])
                        if i==(len(new_metrics)-1):
                            tmpf.write('\n')
                        else:
                            tmpf.write(';')
                f.close()
                tmpf.close()
                df1 = pd.read_csv(outname, sep=';',skiprows=1).dropna()
                df1 = df1[df1.total_score != 'total_score'].drop(columns=['SCORE:'])
                df1['seed']=str(seedid)
                df1['array']=str(array)
                df1.shape
                df = pd.concat([df,df1])
                os.remove(outname)
else:
    raise ValueError("The seed type must be an helix ('H') or a sheet ('S').")

# Formating and writing the output file

df=df[['total_score','buried_unsat_hbonds','buried_unsat_hbonds2','sasa','sc','hbonds','ddG','graft_RMSD','array','seed','description']]

for column in list(df):
    if df[column].dtypes=='float64':
        df[column] = df[column].map('{:.3f}'.format)

df.to_csv('all_score.sc',sep='\t',index=False)
