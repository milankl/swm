## EKE SPEC PLOT
from __future__ import print_function

path1 = '/network/aopp/cirrus/pred/kloewer/swm_back_ronew/'
path2 = '/network/aopp/cirrus/pred/kloewer/swm_bf_cntrl/data/'
outpath = '/network/home/aopp/kloewer/swm/paperplot/'

import os; os.chdir(path2) # change working directory
import numpy as np
from scipy import sparse
import time as tictoc
from netCDF4 import Dataset
import glob
import matplotlib.pyplot as plt
from cmocean import cm

plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'

# OPTIONS
runfolder = [0,6,0,3,8]
print('Compare ek spec plots from run ' + str(runfolder))

runpath1 = path2+'run%04i' % runfolder[0]
D1 = np.load(runpath1+'/analysis/spec_eke.npy').all()
param1 = np.load(runpath1+'/param.npy').all()

runpath2 = path2+'run%04i' % runfolder[1]
D2 = np.load(runpath2+'/analysis/spec_eke.npy').all()
param2 = np.load(runpath2+'/param.npy').all()

runpath3 = path1+'run%04i' % runfolder[2]
D3 = np.load(runpath3+'/analysis/spec_eke.npy').all()
param3 = np.load(runpath3+'/param.npy').all()

runpath4 = path1+'run%04i' % runfolder[3]
D4 = np.load(runpath4+'/analysis/spec_eke.npy').all()
param4 = np.load(runpath4+'/param.npy').all()

runpath5 = path1+'run%04i' % runfolder[4]
D5 = np.load(runpath5+'/analysis/spec_eke.npy').all()
param5 = np.load(runpath5+'/param.npy').all()

Ro_max = param1['c_phase']/(param1['f_0'] - param1['beta']*param1['Ly']/2.)
Ro_min = param1['c_phase']/(param1['f_0'] + param1['beta']*param1['Ly']/2.)

## PLOTTING

kf = 1e3

fig,ax = plt.subplots(1,1,figsize=(7,5))

ax.loglog(D1['k']*kf,D1['p']/(2*np.pi),'C0',label=r'Low resolution $\Delta x = $30km',lw=2)
ax.loglog(D2['k']*kf,D2['p']/(2*np.pi),'C2',label=r'High resolution $\Delta x = $7.5km',lw=2)
ax.loglog(D3['k']*kf,D3['p']/(2*np.pi),'C1',label=r'LR + weak backscatter',ls='--')
ax.loglog(D4['k']*kf,D4['p']/(2*np.pi),'C3',label=r'LR + moderate backscatter',ls='--')
ax.loglog(D5['k']*kf,D5['p']/(2*np.pi),'C5',label=r'LR + strong backscatter',ls='--')

ylim = ax.get_ylim()

ax.loglog(kf/Ro_max*np.ones(2),ylim,'k',lw=0.5)
ax.loglog(kf/Ro_min*np.ones(2),ylim,'k',lw=0.5)

ax.text(1/1900,5e-1,'$L_{Ro}^{max}$',rotation=90,fontsize=12)
ax.text(1/590,5e-1,'$L_{Ro}^{min}$',rotation=90,fontsize=12)

k1 = D1['k'][7]*kf
k2 = D1['k'][-1]*kf

s = 2e-4
ax.loglog([k1,k2],[s*k1**(-3),s*k2**(-3)],'C4',label=r'$K^{-3}$')

xtick = np.array([4000,2000,1000,500,200,100,50,20,10])
xtickm = np.array([3000,900,800,700,600,400,300,90,80,70,60,40,30])
ax.set_xticks(1./xtick)
ax.set_xticks(1./xtickm,minor=True)
ax.set_xticklabels([],minor=True)
ax.set_xticklabels(xtick)

ax.set_xlim(D1['k'][1]*kf,1./xtick[-1])
ax.set_ylim(*ylim)

ax.set_xticklabels(ax.xaxis.get_majorticklabels(), rotation=90)

ax.set_ylabel('EKE [m$^3$s$^{-2}$]')
ax.set_xlabel('wavelength [km]')

ax.legend(loc=1)
ax.set_title('Eddy kinetic energy spectrum',loc='left')

plt.tight_layout()
plt.savefig(outpath+'plots/spec_eke.eps')
plt.close(fig)