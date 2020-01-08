# code to read simulated data from a run
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import os


# scenario
scn = 'reference'
# setting directory
os.chdir('/users/loulou/rise/output/simulation')
# variable for which plots are required by hlth
varnames = ['id','age','wealth','work','claim','income','insurance','ame',
            'cons', 'dead','health','base_health','oop','medexp','value','tr']

# load data
df = pd.read_csv('simulated_'+scn+'.csv',header=None,names=varnames,sep=',',index_col=False)
df = df[(df.age<=84) & (df.age>=35)]
dfw = df[(df.age>=35) & (df.age<=69)]
# recode work
dfw['work'] = dfw['work'].replace(2,0)
df['insurance'] = df['insurance'].replace(2,0)
df['claim'] = df['claim'].replace(1,0)
df['claim'] = df['claim'].replace(2,1)
df['insurance'] = df['insurance'].replace(3,0)
df['wealth'] = df['wealth']/1000.0
df['medexp'] = df['medexp']/1000.0
df['vgood'] = df['base_health'].replace(2,0).replace(3,0).replace(4,1)
df['poor'] = df['base_health'].replace(2,1).replace(3,0).replace(4,0)
df['dead'] = df['dead'].replace(1,0).replace(2,1)

# compute means by age and health
dfh = dfw.groupby(['age','base_health']).mean()
dfa = df.groupby(['age']).mean()

# fit data
varnames = ['data','se']
fit = pd.read_csv('../../params/input/moments.csv',names=varnames)
fit['low'] = fit['data'] - 1.96*fit['se']
fit['high'] = fit['data'] + 1.96*fit['se']
nagesw = 69 - 35 + 1
nages = 84 - 35 + 1
# split for work
work = fit[0:3*nagesw]
pt = 3*nagesw
plt.figure(1,figsize=(9,5),dpi=300)
colors = ['b','r','g']
h = 2
sim = dfh['work'].xs(h,level='base_health')
age = sim.index.values
ww = 0
data = work[ww:ww + nagesw]
ww += nagesw
plt.plot(age,sim,label='simulated', color=colors[h-2])
plt.plot(age,data['data'],label='PSID', linestyle='dashed',color=colors[h-2])
plt.plot(age,data['low'],label='PSID +- 2 s.e.s.', linestyle='-.', color=colors[h-2])
plt.plot(age,data['high'],label='', linestyle='-.', color=colors[h-2])
for h in range(3,5) :
    sim = dfh['work'].xs(h,level='base_health')
    age = sim.index.values
    data = work[ww:ww + nagesw]
    ww += nagesw
    plt.plot(age,sim,label='', color=colors[h-2])
    plt.plot(age,data['data'],label='', linestyle='dashed',color=colors[h-2])
    plt.plot(age,data['low'],label='', linestyle='-.', color=colors[h-2])
    plt.plot(age,data['high'],label='', linestyle='-.', color=colors[h-2])
plt.xlim([35,69])
plt.ylim([0,1])
plt.xticks(range(35,75,5))
plt.legend(loc=3)
plt.xlabel('age')
plt.ylabel('fraction')
plt.title('fraction working')
plt.savefig('../../figures/fit_work.png')
plt.close()

# graph for medexp
medexp = fit[pt:pt + nages]
pt += nages
plt.figure(3,figsize=(9,5),dpi=300)
sim = dfa['medexp']
age = sim.index.values
data = medexp

plt.plot(age,sim,label='simulated', color='b')
plt.plot(age,data['data']/1000.0,label='MEPS', linestyle='dashed',color='b')
plt.plot(age,data['low']/1000.0,label='MEPS +- 2 s.e.s.', linestyle='-.', color='b')
plt.plot(age,data['high']/1000.0,label='', linestyle='-.', color='b')
plt.xlim([35,85])
plt.xticks(range(35,90,5))
plt.ylim([0,30])
plt.legend(loc=2)
plt.xlabel('age')
plt.ylabel('Average medical exp. (000)')
plt.title('Average Medical Expenditures')
plt.savefig('../../figures/fit_medexp.png')
plt.close()



# graph for dead
dead = fit[pt:pt+nages]
pt += nages
plt.figure(4,figsize=(9,5),dpi=300)
sim = dfa['dead']
age = sim.index.values
data = dead
plt.plot(age,sim,label='simulated', color='b')
plt.plot(age,data['data'],label='SSA', linestyle='dashed',color='b')
plt.plot(age,data['low'],label='SSA +- 2 s.e.s.', linestyle='-.', color='b')
plt.plot(age,data['high'],label='', linestyle='-.', color='b')
plt.xlim([35,85])
plt.ylim([0,0.2])
plt.xticks(range(35,90,5))
plt.legend(loc=2)
plt.xlabel('age')
plt.ylabel('fraction')
plt.title('Mortality Rate')
plt.savefig('../../figures/fit_mortality.png')
plt.close()



# graph for poor
poor = fit[pt:pt+nages]
pt += nages
plt.figure(4,figsize=(9,5),dpi=300)
sim = dfa['poor']
age = sim.index.values
data = poor
plt.plot(age,sim,label='simulated', color='b')
plt.plot(age,data['data'],label='PSID', linestyle='dashed',color='b')
plt.plot(age,data['low'],label='PSID +- 2 s.e.s.', linestyle='-.', color='b')
plt.plot(age,data['high'],label='', linestyle='-.', color='b')
plt.xlim([35,85])
plt.ylim([0,0.5])
plt.xticks(range(35,90,5))
plt.legend(loc=2)
plt.xlabel('age')
plt.ylabel('fraction')
plt.title('Fraction Poor/Fair Health')
plt.savefig('../../figures/fit_poor.png')
plt.close()


# graph for vgood
vgood = fit[pt:pt+nages]
pt += nages
plt.figure(4,figsize=(9,5),dpi=300)
sim = dfa['vgood']
age = sim.index.values
data = vgood
plt.plot(age,sim,label='simulated', color='b')
plt.plot(age,data['data'],label='PSID', linestyle='dashed',color='b')
plt.plot(age,data['low'],label='PSID +- 2 s.e.s.', linestyle='-.', color='b')
plt.plot(age,data['high'],label='', linestyle='-.', color='b')
plt.xlim([35,85])
plt.ylim([0,1])
plt.xticks(range(35,90,5))
plt.legend(loc=2)
plt.xlabel('age')
plt.ylabel('fraction')
plt.title('Fraction Very Good/Excellent Health')
plt.savefig('../../figures/fit_vgood.png')
plt.close()

# graph for wealth
wealth = fit[pt:pt+nages]
plt.figure(2,figsize=(9,5),dpi=300)
sim = dfa['wealth']
age = sim.index.values
data = wealth
plt.plot(age,sim,label='simulated', color='b')
plt.plot(age,data['data']/1000.0,label='PSID', linestyle='dashed',color='b')
plt.plot(age,data['low']/1000.0,label='PSID +- 2 s.e.s.', linestyle='-.', color='b')
plt.plot(age,data['high']/1000.0,label='', linestyle='-.', color='b')
plt.xlim([35,85])
plt.ylim([0,750])
plt.xticks(range(35,90,5))
plt.legend(loc=2)
plt.xlabel('age')
plt.ylabel('average wealth (000)')
plt.title('Average Wealth')
plt.savefig('../../figures/fit_wealth.png')
plt.close()

