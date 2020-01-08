# file to run a series of scenarios and collect results
import numpy as np
import os
import pandas as pd
import subprocess
from scipy.optimize import brentq
# change directory where will simulate data
os.chdir('/users/loulou/rise')
# name of scenarios
scenarios = ['full-1965','insure-1965','income-1965','tech-1965','reference']
#scenarios = ['full-1965','insure-1965','income-1965','tech-1965','other-1965','reference']
# var names in simulated datasets
varnames = ['id','age','byear','wealth','work','claim','income','insurance','ame',
            'cons', 'dead','health','base_health','oop','medexp','value','tr']
stats = ['$PDV(c)$','PDV(m)','PDV(y)','$e_{25}$','$e_{50}$','CV(\%)']
labels = ['1965','Insurance','Income','Technology', '2005']

sigma = 3.07
psi = 0.716
beta  = 0.9666
phi_poor = 0.316
phi_good = 0.164
alpha = np.zeros((3,1))
alpha[0] = 0.1611
alpha[1] = alpha[0] + np.exp(-5.309)
alpha[2] = alpha[1] + np.exp(-6.7602)
leisure_endow = 3102.23
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
rate = 0.04

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
    if (work>0.0):
        leisure = leisure - hours_worked*work
    x = (cons**psi)*(leisure**(1.0 - psi))
    u =  alpha[health-2]*base**(1.0-sigma) + (x**(1.0-sigma))/(1.0 - sigma)
    return u
def bequest(wealth):
    beq = phi * ((wealth + kappa)**(psi*(1.0 - sigma)))/(1.0-sigma)
    return beq
def welfare(cv, cons, wealth, work, prob):
    nages = np.size(cons,0)
    nh = 4
    welf = 0.0
    pstate = np.zeros((nages,nh))
    pstate[0][nh-1] = 1.0
    for a in range(0,nages):
        age = agemin + a
        welf += pstate[a][0]*(beta**a)*bequest(wealth[a,0])
        for h in range(1,nh):
            if (cons[a,h-1]>0.0):
                welf += pstate[a][h] * (beta**a) * utility(cons[a,h-1]*(1.0 - cv),work[a,h-1],h+1)
            if (a < nages - 1):
                for hh in range(1,nh):
                    pstate[a+1][h] += pstate[a][hh]*prob[a][hh][h]
        if (a<nages - 1):
            pstate[a+1][0] = np.sum(pstate[a+1][1:nh-1])
    return welf
def solve(cv, cons, wealth, work, prob,valfix):
    return valfix - welfare(cv, cons, wealth, work, prob)
def edpv(var, prob):
    nages = np.size(var,0)
    nh = 4
    pstate = np.zeros((nages,nh))
    pstate[0][nh-1] = 1.0
    value = 0.0
    annuity = 0.0
    for a in range(0,nages):
        for h in range(1,nh):
            annuity += (1.0 - pstate[a][0]) * pstate[a][h] * (1 + rate) ** (-a)
            if (var[a,h-1]>0.0): 
                value += (1.0-pstate[a][0])*pstate[a][h]*var[a,h-1]*(1+rate)**(-a)
            if (a < nages - 1):
                for hh in range(1,nh):
                    pstate[a+1][h] += pstate[a][hh]*prob[a][hh][h]
        if (a<nages - 1):
            pstate[a+1][0] = np.sum(pstate[a+1][1:nh-1])		
    return value/annuity

# number of scenarios
nscn   = len(scenarios)
nstats = 6
results = np.zeros((nscn,nstats))
# flag if want to generate sim data for scenarios
irun = 0
s = 0
trans = np.zeros((nages,nh,nh))
pstate = np.zeros((nages,nh))
cons = np.zeros((nages,nh-1))
work = np.zeros((nages,nh-1))
medexp = np.zeros((nages,nh-1))
wealth = np.zeros((nages,nh-1))
income = np.zeros((nages,nh-1))
for scn in scenarios :
    # load data
    df = pd.read_csv('output/simulation/simulated_cohort_'+scn+'.csv',header=None,names=varnames,sep=',',index_col=False)
    df['work'] = df['work'].replace(2,0)
    # get e25
    sx = df['dead'].groupby(df['age']).size().astype(float)
    sx = sx / sx[sx.index==25].values
    results[s,3] = sx.sum()
    sx50 = sx[sx.index>=50]/sx[sx.index==50].values
    results[s,4] = sx50.sum()

    # compute transition matrix across states
    dums = pd.get_dummies(df['health'],prefix='health')
    df['dead'] = dums['health_1']
    df['poor'] = dums['health_2']
    df['good'] = dums['health_3']
    df['vgood'] = dums['health_4']
    # we do not take into account earnings risk
    for a in range(0,nages):
        age = agemin + a
        # health transition matrix
        for h in range(0,nh) :
            hh = 0
            if (h==0):
                trans[a][h][h] = 1.0
            else :
                for dest in ['dead','poor','good','vgood']:
                    pop = df[(df['base_health']==(h+1)) & (df['age']==age)]
                    if (len(pop.index)!=0):
                        trans[a][h][hh] = pop[dest].mean()
                    hh += 1
            pop = df[(df['age']==age) & (df['base_health']==(h+1))]
            if (len(pop.index)!=0):
                cons[a,h-1] = pop['cons'].mean()
                work[a,h-1] = pop['work'].mean()
                wealth[a,h-1] = pop['wealth'].mean()
                medexp[a,h-1] = pop['medexp'].mean()
                income[a,h-1] = pop['income'].mean()
            else :
                cons[a,h-1] = xmin
                work[a,h-1] = 0.0
                wealth[a,h-1] = 0.0
                medexp[a,h-1] = 0.0
                income[a,h-1] = xmin
            h +=1
    # edpv
    results[s,0] = edpv(cons,trans)
    results[s,1] = edpv(medexp,trans)
    results[s,2] = edpv(income,trans)
    # if scenario is not full-1965: find CV
    if (scn=='full-1965'):
        value1965 = welfare(0.0,cons,wealth,work,trans)
    else :
        cv = brentq(solve, 0.0, 0.8, args=(cons,wealth,work,trans,value1965))
        print('cv = ', cv)
        results[s,5] = cv
    s = s+1

print(results)


results[:,0] = results[:,0]/1e3
results[:,1] = results[:,1]/1e3
results[:,2] = results[:,2]/1e3


# write table to file
table = pd.DataFrame(results,columns=stats, index=scenarios)
tex = open('tables/table-counterfactual.tex', "w")
tex.write('\\begin{tabular}{lrrrrrr} \n')
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
                buff = buff + str(round(results[j,i],3))
        if (i < (nstats-1)) :
            buff = buff + ' & '
    buff = buff + ' \\\ \n'
    tex.write(buff)
tex.write('\hline\hline \n')	
tex.write('\\end{tabular} \n')
tex.close()

#life.get('mx',1991,65)
