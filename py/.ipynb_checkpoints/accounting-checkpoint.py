# file to run a series of scenarios and collect results
import numpy as np
import os
import pandas as pd
import subprocess

# change directory where will simulate data
os.chdir('/users/loulou/rise')
# name of scenarios
scenarios = ['full-1965','insure-1965','income-1965','tech-1965','other-1965','insure-tech-1965','income-tech-1965','income-insure-1965','reference']
scenarios = ['reference']

#scenarios = ['full-1965','insure-1965','income-1965','tech-1965','other-1965','reference']
# var names in simulated datasets
varnames = ['id','age','wealth','work','claim','income','insurance','ame',
            'cons', 'dead','health','base_health','oop','medexp','value','tr']
stats = ['$\overline{m}$','$\overline{oop}$','$e_{25}$','$e_{50}$','CV(\%)']
labels = ['1965','Insurance','Income','Technology','Other','Insurance+Technology','Income+Technology','Income+Insurance','2005']
labshares = ['Insurance','Income','Technology','Other','Insurance+Technology','Income+Technology','Income+Insurance']

sigma = 3.07
psi = 0.715
beta  = 0.9665
phi_poor = 0.316
phi_good = 0.168
alpha = np.zeros((3,1))
alpha[0] = 0.161
alpha[1] = alpha[0] + np.exp(-5.31)
alpha[2] = alpha[1] + np.exp(-6.76)
leisure_endow = 3102.22
hours_worked = 2000.0
phi_poor = (leisure_endow - hours_worked)*phi_poor
phi_good = (leisure_endow - hours_worked)*phi_good
nages = 96
agemin = 25
agemax = 120
nh = 4
phi = 2.529
kappa = 5.0e5
xmin = 13.756e3

# utility function (to compute compensating variations)
def utility(cons, work, health):
    leisure = leisure_endow
    base = (xmin**psi)*(leisure_endow**(1.0 - psi))
    if (health==2):
        leisure = leisure - phi_poor
        base = (xmin**psi)*((leisure_endow - phi_poor)**(1.0 - psi))

    if (health==3):
        leisure = leisure - phi_good
        base = (xmin**psi)*((leisure_endow - phi_good)**(1.0 - psi))
    if (work==1):
        leisure = leisure - hours_worked
    x = (cons**psi)*(leisure**(1.0 - psi))
    u =  alpha[health-2]*base**(1.0-sigma) + (x**(1.0-sigma))/(1.0 - sigma)
    return u
def bequest(wealth):
	beq = phi * ((wealth + kappa)**(psi*(1.0 - sigma)))/(1.0-sigma)
	return beq

def welfare (cons, wealth, prob, cv):
	nages = len(cons)
	nh = 4
	welf = 0.0
	pstate = np.zeros((nages,nh))
	pstate[0][nh-1] = 1.0
	for a in range(0,nages):
		age = agemin + a
		if (age < 65):
			work = 1
		else: 
			work = 0
		welf += pstate[a][0]*(beta**a)*bequest(wealth[a])	
		for h in range(1,nh):
			if (cons[a]>0.0):
				welf += pstate[a][h] * (beta**a) * utility(cons[a]*(1.0 - cv),work,h+1)
			if (a < nages - 1):
				for hh in range(1,nh):
					pstate[a+1][h] += pstate[a][hh]*prob[a][hh][h] 
		if (a<nages - 1):
			pstate[a+1][0] = np.sum(pstate[a+1][1:nh-1])		
		
	return welf


# number of scenarios
nscn   = len(scenarios)
nstats = 5
results = np.zeros((nscn,nstats))
# flag if want to generate sim data for scenarios
irun = 0
s = 0
trans = np.zeros((nscn,nages,nh,nh))
pstate = np.zeros((nscn,nages,nh))
cons = np.zeros((nscn,nages))
wealth = np.zeros((nscn,nages))

for scn in scenarios :
	# running scenario
	if (irun==1) :
		# Calling executable 
		fortran = 'mpirun  --hostfile mpi_hosts.txt ./generate '+scn
		ps = subprocess.Popen(fortran,shell=True,stdin=subprocess.PIPE)
		ps.communicate() #now wait
	# load data
	df = pd.read_csv('../output/simulation/simulated_'+scn+'.csv',header=None,names=varnames,sep=',',index_col=False)
	df['work'] = df['work'].replace(2,0)
	# get spending measures
	results[s,0] = df['medexp'].mean()
	results[s,1] = df['oop'].mean()
	# get e25
	sx = df['dead'].groupby(df['age']).size().astype(float)
	sx = sx / sx[sx.index==25].values
	results[s,2] = sx.sum()
	sx50 = sx[sx.index>=50]/sx[sx.index==50].values
	results[s,3] = sx50.sum()

	# compute transition matrix across states
	dums = pd.get_dummies(df['health'],prefix='health')
	df['dead'] = dums['health_1']
	df['poor'] = dums['health_2']
	df['good'] = dums['health_3']
	df['vgood'] = dums['health_4']	
	for a in range(0,nages):
		age = agemin + a
		# health transition matrix
		for h in range(0,nh) :
			hh = 0
			if (h==0):
				trans[s][a][h][h] = 1.0
			else :	
				for dest in ['dead','poor','good','vgood'] :
					pop = df[(df['base_health']==(h+1)) & (df['age']==age)]
					if (len(pop.index)!=0):
						trans[s][a][h][hh] = pop[dest].mean()						
					hh += 1	
			h +=1

		# consumption path
		pop = df[(df['age']==age)]
		if (len(pop.index)!=0):
			cons[s,a] = pop['cons'].mean()
		else :
			cons[s,a] = 0.0	
		# wealth path
		if (len(pop.index)!=0):
			wealth[s,a] = pop['wealth'].mean()
		else :
			wealth[s,a] = 0.0	
	# if scenario is not full-1965: find CV
	if (scn=='full-1965'):
		value1965 = welfare(cons[s][:],wealth[s][:],trans[s][:][:][:],0.0)
		results[s,4] = 0.0
	else :
		cv = 0.0
		value = 1.0e5
		while (value > value1965):
			value = welfare(cons[s][:],wealth[s][:],trans[s][:][:],cv)
			cv += 0.001
		results[s,4] = cv	

	s = s+1


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
tex = open('../tables/table-counterfactual.tex', "w")
tex.write('\\begin{tabular}{lrrrrr} \n')
tex.write('\hline\hline \n')
buff = 'Outcome & '
for i in range(0,nstats) :
	buff = buff + stats[i] 
	if (i < nstats-1) :
		buff = buff + ' & '
buff = buff + ' \\\ \n'	
tex.write(buff)
tex.write('\hline \n')
for j in range(0,nscn) :
	buff = labels[j] + '	& '
	for i in range(0,nstats) :
		if (i!=nstats-1):
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

#life.get('mx',1991,65)
