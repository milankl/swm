## EKE SPEC PLOT
from __future__ import print_function
path = '/home/mkloewer/python/swm/'
import os; os.chdir(path) # change working directory
import numpy as np
from scipy import sparse
import time as tictoc
from netCDF4 import Dataset
import glob
import matplotlib.pyplot as plt
from cmocean import cm

# OPTIONS
runfolder = [3,10]
print('Compare ek spec plots from run ' + str(runfolder))

runpath1 = path+'data/run%04i' % runfolder[0]
D1 = np.load(runpath1+'/analysis/spec_eke.npy').all()
param1 = np.load(runpath1+'/param.npy').all()

runpath1 = path+'data/run%04i' % runfolder[0]
D11 = np.load(runpath1+'/analysis/spec_eke_b.npy').all()
param1 = np.load(runpath1+'/param.npy').all()

runpath2 = path+'data/run%04i' % runfolder[1]
D2 = np.load(runpath2+'/analysis/spec_eke.npy').all()
param2 = np.load(runpath2+'/param.npy').all()

runpath2 = path+'data/run%04i' % runfolder[1]
D21 = np.load(runpath2+'/analysis/spec_eke_b.npy').all()
param2 = np.load(runpath2+'/param.npy').all()

## PLOTTING

kf = 1e3

fig,ax = plt.subplots(1,1,figsize=(8,6))

ax.loglog(D1['k']*kf,D1['p'],'C0',label=r'Low resolution $\Delta x = 30$km',lw=1.5)
ax.loglog(D11['k']*kf,D11['p'],'C0--',label=r'with boundaries',lw=1.5)
ax.loglog(D2['k']*kf,D2['p'],'C2',label=r'High resolution $\Delta x = 7.5$km',lw=1.5)
ax.loglog(D21['k']*kf,D21['p'],'C2--',label=r'with boundaries',lw=1.5)

k1 = D1['k'][5]*kf
k2 = D1['k'][-1]*kf

s = 2e-3
ax.loglog([k1,k2],[s*k1**(-3),s*k2**(-3)],'C3',label=r'$K^{-3}$')

xtick = np.array([4000,2000,1000,500,200,100,50,20,10])
ax.set_xticks(1./xtick)
ax.set_xticklabels(xtick)

ax.set_xlim(D1['k'][1]*kf,1./xtick[-1])

ax.set_xticklabels(ax.xaxis.get_majorticklabels(), rotation=90)

ax.set_ylabel('EKE [m$^3$s$^{-2}$]')
ax.set_xlabel('wavelength [km]')

ax.legend(loc=1)
ax.set_title('Eddy kinetic energy spectrum')

plt.tight_layout()
plt.savefig(path+'/compare/spec_eke.png')
plt.close(fig)