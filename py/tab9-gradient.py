# file to run a series of scenarios and collect results
import numpy as np
import os
import pandas as pd
import subprocess
from scipy.interpolate import interp1d

# change directory where will simulate data
os.chdir('/home/bepp/michaudp/rise/runtime')
# name of scenarios
# var names in simulated datasets
varnames = ['id','age','byear','wealth','work','claim','income','insurance','ame',
            'cons', 'dead','health','base_health','oop','medexp','value','tr']

# scenario
# setting directory
os.chdir('/home/bepp/michaudp/rise/output/simulation/')

df = pd.read_csv('simulated_pop_reference.csv',header=None,names=varnames,sep=',',index_col=False)
df['work'] = df['work'].replace(2,0) 
df['year'] = df['byear'] + df['age']
df['vgood'] = df['base_health'].replace(2,0).replace(3,0).replace(4,1)
df['poor'] = df['base_health'].replace(2,1).replace(3,0).replace(4,0)
df['good'] = df['base_health'].replace(3,1).replace(2,0).replace(4,0)
df['dead'] = df['dead'].replace(1,0).replace(2,1)
df['qw'] = pd.qcut(df['wealth'],bins=[0.0,100e3,350e3,500e3,1e6],labels=False)
print(df['qw'].value_counts())
df = df[(df['age']>=35) & (df['age']<=85)]
# create first age distribution
ages = [x for x in range(35,85)]
table = pd.DataFrame(index=ages,data=None,columns=['0','1','2','3'])
# get for first quartile 
for q in range(4):
    qdf = df[df['qw']==q]
    pop = qdf['poor'].groupby([qdf['age'],qdf['year']]).mean()
    pop = pop.unstack()
    pop = pop[2005]
    pop = pop[pop.notnull()]
    g = interp1d(pop.index,pop,bounds_error=False,fill_value=pop[-1:])
    poor = g(ages)
    table[str(q)] = poor
table = table[(table.index>=55) & (table.index<=65)]
print(table.mean())







