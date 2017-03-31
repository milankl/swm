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
D11 = np.load(runpath1+'/analysis/spec_eke_b.npy').all()
param1 = np.load(runpath1+'/param.npy').all()

runpath2 = path+'data/run%04i' % runfolder[1]
D2 = np.load(runpath2+'/analysis/spec_eke.npy').all()
D21 = np.load(runpath2+'/analysis/spec_eke_b.npy').all()
param2 = np.load(runpath2+'/param.npy').all()

runpath3 = path+'stoch/data/run%04i' % 0
D3 = np.load(runpath3+'/analysis/spec_eke.npy').all()
D31 = np.load(runpath3+'/analysis/spec_eke_b.npy').all()
param3 = np.load(runpath3+'/param.npy').all()

Ro_max = param1['c_phase']/(param1['f_0'] - param1['beta']*param1['Ly']/2.)
Ro_min = param1['c_phase']/(param1['f_0'] + param1['beta']*param1['Ly']/2.)

## PLOTTING

kf = 1e3

fig,ax = plt.subplots(1,1,figsize=(8,6))

ax.loglog(D1['k']*kf,D1['p'],'C0',label=r'Low resolution $\Delta x = 30$km',lw=1.5)
ax.loglog(D11['k']*kf,D11['p'],'C0--',lw=1.5)
ax.loglog(D2['k']*kf,D2['p'],'C2',label=r'High resolution $\Delta x = 7.5$km',lw=1.5)
ax.loglog(D21['k']*kf,D21['p'],'C2--',lw=1.5)
ax.loglog(D3['k']*kf,D3['p'],'C1',label=r'Low resolution + backscatter',lw=1.5)
ax.loglog(D31['k']*kf,D31['p'],'C1--',lw=1.5)

ylim = ax.get_ylim()

ax.loglog(kf/Ro_max*np.ones(2),ylim,'k',lw=0.5)
ax.loglog(kf/Ro_min*np.ones(2),ylim,'k',lw=0.5)

ax.text(1/1900,5,'$Ro_{max}$',rotation=90,fontsize=12)
ax.text(1/590,5,'$Ro_{min}$',rotation=90,fontsize=12)

ax.loglog(0,0,'0.5',linestyle='dashed',lw=1.5,label='with boundary conditions')

k1 = D1['k'][7]*kf
k2 = D1['k'][-1]*kf

s = 2e-3
ax.loglog([k1,k2],[s*k1**(-3),s*k2**(-3)],'C3',label=r'$K^{-3}$')

xtick = np.array([4000,2000,1000,500,200,100,50,20,10])
ax.set_xticks(1./xtick)
ax.set_xticklabels(xtick)

ax.set_xlim(D1['k'][1]*kf,1./xtick[-1])
ax.set_ylim(*ylim)

ax.set_xticklabels(ax.xaxis.get_majorticklabels(), rotation=90)

ax.set_ylabel('EKE [m$^3$s$^{-2}$]')
ax.set_xlabel('wavelength [km]')

ax.legend(loc=1)
ax.set_title('Eddy kinetic energy spectrum')

plt.tight_layout()
plt.savefig(path+'/compare/spec_eke.png',dpi=300)
plt.close(fig)