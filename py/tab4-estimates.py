# code to produce a table of estimation results
from matplotlib import pyplot as plt 
import pandas as pd
import numpy as np 
import os

# scenarios 

scn = ['reference','nobequest','noleisure','nogrowth','lowk','highsig','highbeta']
scnnames = ['Baseline','No Bequest','No Penalty','No Progress','$K=250e3$','High $\\sigma$','High $\\beta$']
par = ['sigma','psi','phi_good','phi_poor','beta','alpha_poor','alpha_good','alpha_vgood','kapa_base','leisure_endow','prog_theta']
parnames = ['$\sigma$','$\gamma$','$\lambda_1$','$\lambda_2$','$\\beta$','$\\varphi_1$','$\\mu_2$','$\\mu_3$','$\\xi$','$L$','$\\kappa$']
stats = ['OID','degrees','pvalue']
statsname = ['Criterion','D.F.','P-value']
pp = ['est.','st.err']
# setting directory
os.chdir('/home/bepp/michaudp/rise/tables')
n = len(par)+len(stats)
c = len(scn)

pe = pd.DataFrame([], index = par, columns=scn)
se = pd.DataFrame([], index = par, columns=scn)
st = pd.DataFrame([], index = stats, columns=scn)

for sc in scn:
	f = open('../output/estimation/control_estimate_'+sc+'.log')
	lines = f.readlines()
	# results in last X lines
	n = len(lines)
	results = lines[-43:-3]
	for line in results:
		words = line.split()
		if any(words[0] in s for s in par):
			pe[sc].loc[words[0]] = words[1]
			if words[2]=='(fixed)':
				se[sc].loc[words[0]] = '-'
			else :
				se[sc].loc[words[0]] = '('+words[2]+')'
		if any(words[0] in s for s in stats):
			st[sc].loc[words[0]] = words[len(words)-1]

tex = open('table-estimates.tex', "w")
col = 'l'
for r in range(0,len(scn)):
	col = col + 'r' 	
tex.write('\\begin{tabular}{'+col+'} \n')
tex.write('\hline\hline \n')
buff = 'Parameter & '
ii = 0
for i in scn :
	buff = buff + scnnames[ii]
	if (ii < len(scn)-1):
		buff = buff + ' & '
		ii +=1	
buff = buff + ' \\\ \n'	
tex.write(buff)
tex.write('\hline \n')
jj = 0
for j in par :
	buff = parnames[jj] + ' & '
	ii = 0
	for i in scn:
		if (pd.notnull(pe[i].loc[j])):
			buff = buff + str(pe[i].loc[j])	
		else :
			buff = buff + '-'	
			# check 
		if (ii < len(scn)-1): buff = buff + ' & '
		ii += 1
	buff = buff + ' \\\ \n'	
	tex.write(buff)
	buff = ' & '
	ii = 0
	for i in scn:
		if (pd.notnull(pe[i].loc[j])):
			buff = buff + str(se[i].loc[j])	
		else :
			buff = buff + '-'	
			# check 
		if (ii < len(scn)-1): buff = buff + ' & '
		ii += 1
	buff = buff + ' \\\ \n'	
	tex.write(buff)
	jj += 1
tex.write('\hline \n')
jj = 0	
for j in stats :
	buff = statsname[jj] + ' & '
	ii = 0
	for i in scn:
		if (pd.notnull(st[i].loc[j])):
			buff = buff + str(st[i].loc[j])	
		else :
			buff = buff + '-'	
			# check 
		if (ii < len(scn)-1): buff = buff + ' & '
		ii += 1
	buff = buff + ' \\\ \n'	
	tex.write(buff)	
	jj +=1
tex.write('\hline\hline \n')	
tex.write('\\end{tabular} \n')
tex.close()









