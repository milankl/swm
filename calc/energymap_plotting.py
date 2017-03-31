## PLOT ENERGY RESERVOIR MAPS
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

##
"""
for v in ['mke','eke','epe','mpe']:
    D1[v] = np.log10(D1[v])
    D2[v] = np.log10(D2[v])
"""
## PLOTTING   
levs=[]
for v in ['mke','eke','mpe','epe']:
    maxs = np.max((np.max(D1[v]),np.max(D2[v])))
    mins = np.min((np.min(D1[v]),np.min(D2[v])))
    levs.append(np.linspace(mins,maxs*0.8,64))

fig,axs = plt.subplots(2,4,figsize=(14,8),sharex=True,sharey=True)
plt.tight_layout(rect=[0,.08,1,1])
n = axs.shape[1]
caxs = [0,]*n
for i in range(n):
    pos = axs[-1,i].get_position()
    caxs[i] = fig.add_axes([pos.x0,0.05,pos.width,0.03])

qaxs = np.empty_like(axs)

qaxs[0,0] = axs[0,0].contourf(param1['x_T']/1e3,param1['y_T']/1e3,h2mat(D1['mke'],param1),levs[0],cmap='viridis',extend='max')
fig.colorbar(qaxs[0,0],cax=caxs[0],orientation='horizontal',ticks=[0,1e5,2e5])

qaxs[0,1] = axs[0,1].contourf(param1['x_T']/1e3,param1['y_T']/1e3,h2mat(D1['eke'],param1),levs[1],cmap='viridis',extend='max')
fig.colorbar(qaxs[0,1],cax=caxs[1],orientation='horizontal',ticks=[0,5e5,1e6])

qaxs[0,2] = axs[0,2].contourf(param1['x_T']/1e3,param1['y_T']/1e3,h2mat(D1['mpe'],param1),levs[2],cmap='viridis',extend='max')
fig.colorbar(qaxs[0,2],cax=caxs[2],orientation='horizontal',ticks=[0,5e3,1e4])

qaxs[0,3] = axs[0,3].contourf(param1['x_T']/1e3,param1['y_T']/1e3,h2mat(D1['epe'],param1),levs[3],cmap='viridis',extend='max')
fig.colorbar(qaxs[0,3],cax=caxs[3],orientation='horizontal',ticks=[0,1e4,4e4])

axs[1,0].contourf(param2['x_T']/1e3,param2['y_T']/1e3,h2mat(D2['mke'],param2),levs[0],cmap='viridis',extend='max')
axs[1,1].contourf(param2['x_T']/1e3,param2['y_T']/1e3,h2mat(D2['eke'],param2),levs[1],cmap='viridis',extend='max')    
axs[1,2].contourf(param2['x_T']/1e3,param2['y_T']/1e3,h2mat(D2['mpe'],param2),levs[2],cmap='viridis',extend='max')    
axs[1,3].contourf(param2['x_T']/1e3,param2['y_T']/1e3,h2mat(D2['epe'],param2),levs[3],cmap='viridis',extend='max')    

axs[0,0].set_title(r'MKE [Jm$^{-2}$]')
axs[0,1].set_title(r'EKE [Jm$^{-2}$]')
axs[0,2].set_title(r'MPE [Jm$^{-2}$]')
axs[0,3].set_title(r'EPE [Jm$^{-2}$]')

axs[0,0].set_yticks([])
axs[1,0].set_yticks([])
axs[1,0].set_xticks([])
axs[1,1].set_xticks([])
axs[1,2].set_xticks([])
axs[1,3].set_xticks([])

axs[0,0].set_ylabel(r'Low resolution, $\Delta x = 30$km')
axs[1,0].set_ylabel(r'High resolution, $\Delta x = 7.5$km')

plt.savefig(path+'compare/Energy_mean.png')
plt.close(fig)
