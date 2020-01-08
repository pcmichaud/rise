# file to run a series of scenarios and collect results
import numpy as np
import os
import pandas as pd
import subprocess
from scipy.interpolate import interp1d

# change directory where will simulate data
os.chdir('/home/bepp/michaudp/rise/runtime')
# name of scenarios
scenarios = ['full-1965','insure-1965','income-1965','tech-1965','other-1965','insure-tech-1965','income-tech-1965','income-insure-1965','reference']
#scenarios = ['reference']
#scenarios = ['full-1965','insure-1965','income-1965','tech-1965','other-1965','reference']
# var names in simulated datasets
varnames = ['id','age','byear','wealth','work','claim','income','insurance','ame',
            'cons', 'dead','health','base_health','oop','medexp','value','tr']
stats = ['$\overline{m}$','$\overline{oop}$','$e_{25}$','$e_{50}$']
labels = ['1965','Insurance','Income','Technology','Other','Insurance+Technology','Income+Technology','Income+Insurance','2005']
labshares = ['Insurance','Income','Technology','Other','Insurance+Technology','Income+Technology','Income+Insurance']

# number of scenarios
nscn   = len(scenarios)
nstats = 4
results = np.zeros((nscn,nstats))
# flag if want to generate sim data for scenarios
irun = 0
s = 0



for scn in scenarios :
    # load data
    df = pd.read_csv('../output/simulation/simulated_pop_'+scn+'.csv',header=None,names=varnames,sep=',',index_col=False)
    df['work'] = df['work'].replace(2,0) 
    df['year'] = df['byear'] + df['age']
    print(df.head())
    # create first age distribution
    ages = [x for x in range(25,100)]
    pop = df['medexp'].groupby([df['age'],df['year']]).count()
    pop = pop.unstack()
    print(pop.columns)
    pop = pop[2005]
    pop = pop[pop.notnull()]
    g = interp1d(pop.index,pop,bounds_error=False,fill_value=pop[-1:])
    freq = g(ages)
    freq = freq/freq.sum()

    # get average medical expenditures
    med = df['medexp'].groupby([df['age'],df['year']]).mean()
    med = med.unstack()
    med = med[2005]
    med = med[med.notnull()]
    f = interp1d(med.index,med,bounds_error=False,fill_value=med[-1:])
    pmed = f(ages)
    results[s,0] = np.sum([m*freq[i] for i,m in enumerate(pmed)])

    # get average out-of-pocket medical expenditures
    med = df['oop'].groupby([df['age'],df['year']]).mean()
    med = med.unstack()
    med = med[2005]
    med = med[med.notnull()]
    f = interp1d(med.index,med,bounds_error=False,fill_value=med[-1:])
    pmed = f(ages)
    results[s,1] = np.sum([m*freq[i] for i,m in enumerate(pmed)])

    df['dead'] = df['dead'].replace(1,0).replace(2,1)
    mx = df['dead'].groupby([df['age'],df['year']]).mean()
    mx = mx.unstack()
    mx = mx[2005]
    mx = mx[mx.notnull()]
    f = interp1d(mx.index,mx,bounds_error=False,fill_value=mx[-1:])
    pmx = f(ages)
    sx = [1.0]
    for i in range(1,len(pmx)):
        sx.append(sx[i-1]*(1-pmx[i-1]))
    results[s,2] = np.sum(sx)
    results[s,3] = np.sum([sx[i]/sx[24] for i in range(24,len(sx))])
    s = s+1

print(results)

# computing share of progress to add to columns for medexp and lifeexp
share = np.zeros((nscn,4))
dtot_med = results[nscn-1,0] - results[0,0]
dtot_oop = results[nscn-1,1] - results[0,1]
dtot_lex = results[nscn-1,2] - results[0,2]
dtot_e50 = results[nscn-1,3] - results[0,3]

for s in range(0,nscn): 
	share[s,0] = (results[s,0] - results[0,0])/dtot_med
	share[s,1] = (results[s,1] - results[0,1])/dtot_oop
	share[s,2] = (results[s,2] - results[0,2])/dtot_lex
	share[s,3] = (results[s,3] - results[0,3])/dtot_e50


# write table to file
table = pd.DataFrame(results,columns=stats, index=scenarios)
tex = open('../tables/table-counterfactual-pop.tex', "w")
tex.write('\\begin{tabular}{lrrrr} \n')
tex.write('\hline\hline \n')
buff = 'Outcome & '
for i in range(0,nstats) :
	buff = buff + stats[i] 
	if (i < nstats) :
		buff = buff + ' & '
buff = buff + ' \\\ \n'	
tex.write(buff)
tex.write('\hline \n')
for j in range(0,nscn) :
	buff = labels[j] + '	& '
	for i in range(0,nstats) :
		if (i!=nstats):
			buff = buff + str(round(results[j,i],1))
		else:
			if (j!=0):
				buff = buff + str(round(results[j,i],2))
		if (i < (nstats-1)) :
			buff = buff + ' & '
	buff = buff + ' \\\ \n'
	tex.write(buff)
	if (any(labels[j] in s for s in labshares)):
		buff = '\\multicolumn{1}{r}{$\Delta$\% Total} & '+ str(round(share[j,0],3)) + ' & ' + str(round(share[j,1],3)) + '  &  ' + str(round(share[j,2],3)) + ' & ' + str(round(share[j,3],3)) + '  &   \\\ \n' 
		tex.write(buff)
tex.write('\hline\hline \n')	
tex.write('\\end{tabular} \n')
tex.close()
