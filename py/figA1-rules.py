    # code to read decision rules
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt 
import pandas as pd
import numpy as np 
import os
from scipy.interpolate import interp1d
# scenario 
scn = 'reference'
# setting directory
os.chdir('/home/bepp/michaudp/rise/output/rules')

# variable for which plots are required by hlth
varnames = ['medexp','cons','value','claim','work']
colnames = ['h','e','z','f','d','ame','w','iwealth','cons','medexp',
            'work','claim','value','cash']

# state chosen
ee = 5
z = 1
dd = 1
f = 3
ame = 12
wmax = 5.0e5

ages = [35, 45, 55, 70, 75, 80, 85]

tex = open('rules_'+scn+'.tex', "w")
tex.write('\\documentclass{beamer} \n')
tex.write('\\title{Decision Rules: Accounting for the Rise of Health Spending and Longevity} \n')
tex.write('\\author{Fonseca, Michaud, Kapteyn and Galama} \n')
tex.write('\\date{\\today} \n')
tex.write('\\begin{document} \n')
tex.write('\\frame{\\titlepage} \n')

# age chosen
s = 1
for age in ages:
    df = pd.read_csv('rules_age'+str(age)+'_'+scn+'.csv',header=None,
                     names = colnames)  
 
    if (age <= 62) :
        d = 1
    elif (age>=70) :
        d = 2
    else :    
        d = dd
    if (age < 70) :
        e = ee
    elif (d == 2) :
        e = 1    
    else :
        e = 1            
    for v in varnames :
        plt.figure(s)  
        tex.write('\\frame{ \n')
        tex.write('\\frametitle{'+v+' at '+str(age)+'} \n')
        for h in range(2,5): 
            cond = (df['h']==h) & (df['e'] == e) & (df['z'] == z) & (df['d'] == d) & (df['f'] == f) & (df['ame'] == ame) & (df['iwealth'] < wmax) & (df['iwealth'] > 1.0)
            data = df[cond]
            var = data[v].values
            wlth = data['iwealth'].values
            fn = interp1d(wlth,var,kind='cubic')
            plt.plot(wlth,fn(wlth),label='health = '+str(h-1))
            plt.xlabel('assets')
            plt.ylabel(v)
            plt.legend(loc=2)
            plt.title('Rules at age '+str(age)+' for : '+v)
            plt.savefig('rules_age'+str(age)+'_'+v+'_'+scn+'.png')
        s = s + 1
        plt.close()
        tex.write('\\begin{figure} \n')
        tex.write('\\includegraphics[scale=0.5]{rules_age'+str(age)+'_'+v+'_'+scn+'.png} \n')
        tex.write('\\caption{rules evaluated at (e,z,d,f,ame) = ('+str(e)+','+str(z)+','+str(d)+','+str(f)+','+str(ame)+').}')        
        tex.write('\\end{figure} \n')
        tex.write('} \n')
        
tex.write('\\end{document} \n')
tex.close()  
#os.system('pdflatex -interaction=batchmode rules_'+scn+'.tex') 
