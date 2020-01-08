# code to read simulated data from a run
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import os
import subprocess
os.chdir('/home/bepp/michaudp/rise/runtime')
# scenario
scenarios = ['copay95','copay75','copay50','copay25']
scnnames = ['$\overline{m}(95)$','$\overline{m}(75)$','$\overline{m}(50)$','$\overline{m}(25)$']
copay = [0.95,0.75,0.5,0.25]
enames = ['e(95)','e(75)','e(50)']
elatex = ['$\epsilon_p(95)$','$\epsilon_p(75)$','$\epsilon_p(50)$']
# setting directory
# variable for which plots are required by hlth
varnames = ['id','age','byear','wealth','work','claim','income','insurance','ame',
            'cons', 'dead','health','base_health','oop','medexp','value','tr']
agenames = ['30-34']
for i in range(35,85,5):
	lab = [str(i)+'-'+str(i+4)]
	agenames = agenames + lab

irun = 0
ss = 0
result =pd.DataFrame(index=agenames,columns=scenarios)

for scn in scenarios :
	# running scenario
	if (irun==1) :
		# Calling executable
		fortran = '/opt/openmpi/intel/bin/mpirun  --hostfile mpi_hosts.txt ./generate '+scn
		ps = subprocess.Popen(fortran,shell=True,stdin=subprocess.PIPE)
		ps.communicate() #now wait
	# load data
	df = pd.read_csv('../output/simulation/simulated_'+scn+'.csv',header=None,names=varnames,sep=',',index_col=False)
	df = df[(df['age']>=30) & (df['age']<85)]
	df['age5'] = pd.cut(df['age'], 11,labels=agenames)
	df['copay'] = copay[ss]
	result[scn] = df.groupby(['age5'])['medexp'].mean()
	if (ss!=0):
		lastscn = scenarios[ss-1]
		c = copay[ss-1]*100
		var = 'e('+str(int(c))+')'
		dp = (copay[ss] - copay[ss-1])/(copay[ss] + copay[ss-1])
		result[var] = (result[scn] - result[lastscn])/(result[scn]+result[lastscn])
		result[var] = result[var]/dp
	ss +=1


tex = open('../tables/table-price-elasticity.tex', "w")
col = 'l'
for r in range(0,len(scenarios)):
	col = col + 'r'
for r in range(0,len(enames)):
	col = col + 'r'
tex.write('\\begin{tabular}{'+col+'} \n')
tex.write('\hline\hline \n')
buff = 'Age & '
ii = 0
for i in scenarios :
	buff = buff + scnnames[ii]
	buff = buff + ' & '
	ii +=1
ii = 0
for i in enames :
	buff = buff + elatex[ii]
	if (ii < len(enames)-1):
		buff = buff + ' & '
		ii +=1
buff = buff + ' \\\ \n'
tex.write(buff)
tex.write('\hline \n')
jj = 0
for j in agenames :
	buff = j + ' & '
	ii = 0
	for i in scenarios:
		buff = buff + str(round(result[i].loc[j],1))
		buff = buff + ' & '
		ii += 1
	ii = 0
	for i in enames:
		buff = buff + str(round(result[i].loc[j],3))
		if (ii < len(enames)-1): buff = buff + ' & '
		ii += 1
	buff = buff + ' \\\ \n'
	tex.write(buff)
tex.write('\hline\hline \n')
tex.write('\\end{tabular} \n')
tex.close()
