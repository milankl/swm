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

nfolders = [0,1,2]

D = []
Db = []
p = []

for r in nfolders:
    runpath = '/home/mkloewer/github/swm/data/run%04i' % r
    D.append(np.load(runpath+'/analysis/spec_eke.npy').all())
    Db.append(np.load(runpath+'/analysis/spec_eke_b.npy').all())
    p.append(np.load(runpath+'/param.npy').all())

Ro_max = param1['c_phase']/(param1['f_0'] - param1['beta']*param1['Ly']/2.)
Ro_min = param1['c_phase']/(param1['f_0'] + param1['beta']*param1['Ly']/2.)

## PLOTTING
kf = 1e3

fig,ax = plt.subplots(1,1,figsize=(8,6))

l1, = ax.loglog(D11['k']*kf,D11['p'],'C0',label=r'Low resolution $\Delta x = 30$km',lw=2.5)
l2, = ax.loglog(D21['k']*kf,D21['p'],'C2',label=r'High resolution $\Delta x = 7.5$km',lw=2.5)

leg1 = ax.legend(handles=[l1,l2],loc=3,ncol=2,title='Control runs')

ll = []
for db,ip,co in zip(Db,p,[3,1,4]):
    pars = (ip['c_D'])
    l, = ax.loglog(db['k']*kf,db['p'],'C'+str(co),label=r'$c_{D}$=%.2e' % pars,lw=1.5)
    ll.append(l)

ylim = ax.get_ylim()

ax.loglog(kf/Ro_max*np.ones(2),ylim,'k',lw=0.5)
ax.loglog(kf/Ro_min*np.ones(2),ylim,'k',lw=0.5)

ax.text(1/1900,25,'$Ro_{max}$',rotation=90,fontsize=12)
ax.text(1/590,25,'$Ro_{min}$',rotation=90,fontsize=12)


k1 = D1['k'][7]*kf
k2 = D1['k'][-1]*kf

s = 2e-3
l, = ax.loglog([k1,k2],[s*k1**(-3),s*k2**(-3)],'C8',label=r'$K^{-3}$')
ll.append(l)

ax.legend(handles=ll,loc=1,title='with backscatter')
ax.add_artist(leg1)


xtick = np.array([4000,2000,1000,500,200,100,50,20,10])
ax.set_xticks(1./xtick)
ax.set_xticklabels(xtick)

ax.set_xlim(D1['k'][1]*kf,1./xtick[-1])
ax.set_ylim(*ylim)

ax.set_xticklabels(ax.xaxis.get_majorticklabels(), rotation=90)

ax.set_ylabel(r'EKE [m$^3$s$^{-2}$]')
ax.set_xlabel('wavelength [km]')

ax.set_title('Eddy kinetic energy spectrum')

plt.tight_layout()
plt.savefig(path+'/compare/spec_eke_bf.png',dpi=150)
plt.close(fig)