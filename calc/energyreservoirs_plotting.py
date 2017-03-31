## PRODUCE ENERGY RESERVOIRS PLOT

from __future__ import print_function
path = '/home/mkloewer/python/swm/'
import os; os.chdir(path) # change working directory
import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt
import time as tictoc
from netCDF4 import Dataset
import glob
from matplotlib.colors import BoundaryNorm,LogNorm
import cmocean

# OPTIONS
runfolder = [3,10]
print('Compare mean plots from run ' + str(runfolder))

## read data

runpath1 = path+'data/run%04i' % runfolder[0]
D1 = np.load(runpath1+'/analysis/mean.npy').all()
param1 = np.load(runpath1+'/param.npy').all()

runpath2 = path+'data/run%04i' % runfolder[1]
D2 = np.load(runpath2+'/analysis/mean.npy').all()
param2 = np.load(runpath2+'/param.npy').all()

# functions
def h2mat(h,param):
    return h.reshape((param['ny'],param['nx']))

def u2mat(u,param):
    return u.reshape((param['ny'],param['nx']-1))

def v2mat(v,param):
    return v.reshape((param['ny']-1,param['nx']))

def q2mat(q,param):
    return q.reshape((param['ny']+1,param['nx']+1))

## PLOTTING   

for v in ['mke','eke','epe','mpe']:
    D1[v] = D1[v].mean()
    D2[v] = D2[v].mean()

me = np.array([D1['mke'],D1['mpe'],D2['mke'],D2['mpe']]) * 1e-3   # mean energies
ee = np.array([D1['eke'],D1['epe'],D2['eke'],D2['epe']]) * 1e-3   # eddy energies

x = np.array([.5,1.5,3.5,4.5])

fig,(ax1,ax2) = plt.subplots(2,1,sharex=True)

ax1.bar(x[::2],me[::2],color='g',alpha=.5,label='kinetic')
ax1.bar(x[1::2],me[1::2],color='b',alpha=.5,label='potential')
ax2.bar(x,ee,color=['g','b','g','b'],alpha=.5)

ax1.legend(loc=5,fontsize=10)

ax1.text(0.25,.9,r'Low resolution, $\Delta x = 30$km',ha='center',transform=ax1.transAxes)
ax2.text(0.25,.9,r'Low resolution, $\Delta x = 30$km',ha='center',transform=ax2.transAxes)
ax1.text(0.75,.9,r'High resolution, $\Delta x = 7.5$km',ha='center',transform=ax1.transAxes)
ax2.text(0.75,.9,r'High resolution, $\Delta x = 7.5$km',ha='center',transform=ax2.transAxes)

ax1.text(x[0]+.4,me[0]*.99,'MKE',color='#444444',ha='center',va='top')
ax1.text(x[1]+.4,me[1],'MPE',color='#444444',ha='center',va='bottom')
ax1.text(x[2]+.4,me[2]*.99,'MKE',color='#444444',ha='center',va='top')
ax1.text(x[3]+.4,me[3],'MPE',color='#444444',ha='center',va='bottom')

ax2.text(x[0]+.4,ee[0]*.99,'EKE',color='#444444',ha='center',va='top')
ax2.text(x[1]+.4,ee[1],'EPE',color='#444444',ha='center',va='bottom')
ax2.text(x[2]+.4,ee[2]*.99,'EKE',color='#444444',ha='center',va='top')
ax2.text(x[3]+.4,ee[3],'EPE',color='#444444',ha='center',va='bottom')

ax1.text(.47,.6,r'$\times$%1.2f' % (me[2]/me[0]),ha='center',transform=ax1.transAxes,color='g')
ax1.arrow(.43,.55,.07,0,transform=ax1.transAxes,fc='g',ec='g',alpha=.7)
ax1.text(.47,.4,r'$\times$%1.2f' % (me[3]/me[1]),ha='center',transform=ax1.transAxes,color='b')
ax1.arrow(.43,.35,.07,0,transform=ax1.transAxes,fc='b',ec='b',alpha=.7)
ax2.text(.47,.6,r'$\times$%1.2f' % (ee[2]/ee[0]),ha='center',transform=ax2.transAxes,color='g')
ax2.arrow(.43,.55,.07,0,transform=ax2.transAxes,fc='g',ec='g',alpha=.7)
ax2.text(.47,.4,r'$\times$%1.2f' % (ee[3]/ee[1]),ha='center',transform=ax2.transAxes,color='b')
ax2.arrow(.43,.35,.07,0,transform=ax2.transAxes,fc='b',ec='b',alpha=.7)


ax1.set_title('Mean energy')
ax2.set_title('Eddy energy')

ax1.set_ylabel(r'kJm$^{-2}$')
ax2.set_ylabel(r'kJm$^{-2}$')

ax1.set_ylim(0,D2['mke']*1.2*1e-3)
ax2.set_xticks([])

plt.tight_layout()
plt.savefig(path+'compare/Energy_reservoirs.png')
plt.close(fig)
#plt.show()