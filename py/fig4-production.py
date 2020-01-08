import numpy as np
from matplotlib import pyplot as plt

# global parameters
nh = 3

# health
delta_p = np.matrix([[0,.6124402,1.266343],[0,-.0206929,-.0400112]])
delta_a = np.array([0,-.0527256,-.1054228])
delta_h = np.matrix([[0,0,0],
    [0,2.527637,3.492382],
    [0,3.148639,6.065002]])
delta_r = np.matrix([[0,-.2603543,-.6744807],
    [0,-.1039402,-.2273313]])
delta_b = np.array([0,-1.062134,-2.744025])
# mortality
psi_a = .0783022
psi_h = np.array([0, -1.147022, -1.886106])
psi_b  = -8.251272

kappa = 0.0071

# progress function
def d(t):
    return np.exp(kappa*(t-2005))

# functions
def i(j,medexp,t):
    return d(t)*(delta_p[0,j]*np.log(1.0+medexp) +
     delta_p[1,j]*np.log(1.0+medexp)**2)
def probh(j,k,medexp,o,s,a,t):
    H = np.zeros(nh)
    for jj in range(nh):
        iv = i(jj,medexp,t)
        H[jj] = delta_a[jj]*a + iv + delta_h[k,jj] + delta_r[0,jj]*o + delta_r[1,jj]*s + delta_b[jj]
        H[jj] = np.exp(H[jj])
    return H[j]/np.sum(H)
def probm(k,a):
    H = psi_a*a + psi_h[k] + psi_b
    return 1.0-np.exp(-np.exp(H))
def probs(k,medexp,o,s,a,t):
    ps = 0.0
    for j in range(nh):
        ph = probh(j,k,medexp,o,s,a,t)
        pm = probm(j,a)
        ps += ph*(1.0-pm)
    return ps
def mcost(medexp,ps):
    n = len(medexp)
    mc = []
    for i in range(n-1):
        dm = (medexp[i+1]-medexp[i])*1e-3
        ds = (ps[i+1] - ps[i])
        if ds>0.0:
            mc.append((dm/ds)/100.0)
        else :
            mc.append(np.nan)
    return mc

# plot
o = 0
s = 0
mgrid = np.linspace(100.0,50e3,100)
labels = ['poor','good','very good']
color=['r','b','g']
plt.figure()
for k in range(nh):
    ps = [probs(k,m,o,s,75,2005) for m in mgrid]
    plt.plot(mgrid*1e-3,ps,color=color[k],label=labels[k]+' - 2005')
    ps = [probs(k,m,o,s,75,1965) for m in mgrid]
    plt.plot(mgrid*1e-3,ps,color=color[k],linestyle='--',label=labels[k]+' - 1965')


plt.xlabel('medical expenditures (thousands)')
plt.ylabel('survival probability')
plt.xlim([0.1,50])
plt.legend()
plt.savefig('../figures/production_change.png')
plt.show()

# plot
o = 0
s = 0
nm = 100
mgrid = np.linspace(100.0,50e3,nm)
labels = ['poor','good','very good']
color=['r','b','g']
plt.figure()
for k in range(nh):
    ps = [probs(k,m,o,s,75,2005) for m in mgrid]
    mc2005 = mcost(mgrid,ps)
    plt.plot(ps[:-1],mc2005,color=color[k],label=labels[k]+' - 2005')
    ps = [probs(k,m,o,s,75,1965) for m in mgrid]
    mc1965 = mcost(mgrid,ps)
    plt.plot(ps[:-1],mc1965,color=color[k],linestyle='--',label=labels[k]+' - 1965')
plt.xlabel('survival probability')
plt.ylabel('marginal KeyboardInterruptcost of increasing survival 1% (000)')
plt.ylim([0.0,250])
plt.legend()
plt.savefig('../figures/cost_change.png')
plt.show()