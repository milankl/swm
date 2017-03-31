## PLOT PV MEAN
from __future__ import print_function
path = '/home/mkloewer/python/swm/'
import os; os.chdir(path) # change working directory
import numpy as np
from scipy import sparse
from scipy.integrate import cumtrapz
import matplotlib.pyplot as plt
import time as tictoc
from netCDF4 import Dataset
import glob
from cmocean import cm

# OPTIONS
runfolder = [3,10]
print('Compare mean plots from run ' + str(runfolder))

## read data

runpath1 = path+'data/run%04i' % runfolder[0]
D1 = np.load(runpath1+'/analysis/PVm_Lm.npy').all()
D11 = np.load(runpath1+'/analysis/mean.npy').all()
param1 = np.load(runpath1+'/param.npy').all()

runpath2 = path+'data/run%04i' % runfolder[1]
D2 = np.load(runpath2+'/analysis/PVm_Lm.npy').all()
D21 = np.load(runpath2+'/analysis/mean.npy').all()
param2 = np.load(runpath2+'/param.npy').all()

#runpath3 = path+'stoch/data/run%04i' % 7
#D3 = np.load(runpath3+'/analysis/PVm_Lm.npy').all()
#param3 = np.load(runpath3+'/param.npy').all()

# import functions
exec(open(path+'swm_param.py').read())
exec(open(path+'swm_operators.py').read())
exec(open(path+'swm_output.py').read())
exec(open(path+'swm_rhs.py').read())
param1['output'] = 0
param2['output'] = 0

## Plotting
fig,axs = plt.subplots(1,2,sharex=True,sharey=True,figsize=(6,4))
fig.tight_layout(rect=[0,.13,1,0.95])
fig.subplots_adjust(wspace=0.05,hspace=0.05)

pos = axs[0].get_position()
pos2 = axs[-1].get_position()
cax = fig.add_axes([pos.x0,0.13,pos2.x1-pos.x0,0.03])

s = 1e7 # scaling factor

levs = np.linspace(np.percentile(D2['PVm'],1)*s,np.percentile(D2['PVm'],99)*s,32)

global param
for i,(D,p) in enumerate(zip([D1,D2],[param1,param2])):
    param = p
    q = axs[i].contourf(p['x_q'],p['y_q'],q2mat(D['PVm'])*1e7,levs,cmap=cm.thermal,extend='both')

cbar = fig.colorbar(q,cax=cax,orientation='horizontal')
cbar.set_label('Potential Vorticity [$m^{-1}s^{-1} \cdot$ 1e-7]')
cbar.set_ticks([1,1.25,1.5,1.75,2])
axs[0].set_xticks([0,p['Lx']])
axs[0].set_yticks([0,p['Ly']])

axs[0].set_xticklabels([0,r'$L_x$'])
axs[0].set_yticklabels([0,r'$L_y$'])

axs[0].set_xlim(0,param['Lx'])
axs[0].set_ylim(0,param['Ly'])

axs[0].set_title(r'Low resolution, $\Delta x = 30$km')
axs[1].set_title(r'High resolution, $\Delta x = 7.5$km')

plt.savefig(path+'compare/pv_mean.png',dpi=300)
plt.close(fig)