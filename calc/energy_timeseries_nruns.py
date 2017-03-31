## PLOT ENERGY AND ENSTROPY TIME SERIES
from __future__ import print_function
path = '/home/mkloewer/python/swm/'
import os; os.chdir(path) # change working directory
import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt
import time as tictoc
from netCDF4 import Dataset
import glob
from cmocean import cm
import matplotlib.patches as patches

# OPTIONS
runfolder = [1,5,6]
print('Compare mean plots from run ' + str(runfolder))

## read data

D = []
params = []

runpath1 = path+'data/run%04i' % runfolder[0]
D.append(np.load(runpath1+'/analysis/mean.npy').all())

runpath2 = path+'data/run%04i' % runfolder[1]
D.append(np.load(runpath2+'/analysis/mean.npy').all())

runpath3 = path+'data/run%04i' % runfolder[2]
D.append(np.load(runpath3+'/analysis/mean.npy').all())

runs = [0,1,2]

for r in runs:
    runpath = '/home/mkloewer/github/swm/data/run%04i' % r
    D.append(np.load(runpath+'/analysis/timeseries.npy').all())
    params.append(np.load(runpath+'/param.npy').all())

## PLOTTING   

fig,ax1 = plt.subplots(1,1,sharex=True,figsize=(8,6))

ys = 365*24*3600    # one year in seconds
ks = 1e-3          # scaling for energy
zs = 1e11           # scaling for enstrophy

# control runs
l1, = ax1.plot(D[0]['time']/ys,D[0]['KEm']*ks,'C0',label=r'Low resolution, $\Delta x = 30$km',lw=2)
l2, = ax1.plot(D[1]['time']/ys,D[1]['KEm']*ks,'C2',label=r'High resolution, $\Delta x = 7.5$km',lw=2)
ax1.plot(D[2]['time']/ys,D[2]['KEm']*ks,'C2',lw=2)

leg1 = ax1.legend(handles=[l1,l2],loc=4,ncol=2,title='Control runs')

ll = []
for i,d,p,ic in zip(runs,D[3:],params,[3,1,4]):
    pars = (p['c_D'])
    l, = ax1.plot(d['time']/ys,d['KEm']*ks,'C'+str(ic),lw=1,label=r'$c_{D}$=%.2e' % pars)
    ll.append(l)


ax1.legend(handles=ll,loc=1,title='with backscatter')
ax1.add_artist(leg1)

ax1.set_title('Kinetic Energy')
ax1.set_ylabel(r'kJm$^{-2}$')

ax1.set_xlim(0,10)
ax1.set_ylim(0,330)

ax1.set_xlabel('years')

plt.tight_layout()
plt.savefig('compare/Energy_timeseries_nruns_bf.png',dpi=200)
#plt.show()
plt.close(fig)