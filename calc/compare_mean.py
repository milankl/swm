## PRODUCE MEAN PLOTS COMPARE FROM TWO RUNS
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
hm_max = np.max(abs(D2['hm']))
hlevs = np.linspace(-hm_max,hm_max,64)

um_max = np.max(abs(D2['um']))
ulevs = np.linspace(-um_max,um_max,64)

vm_max = np.max(abs(D2['vm']))
vlevs = np.linspace(-vm_max,vm_max,64)

fig,((ax11,ax12,ax13),(ax21,ax22,ax23)) = plt.subplots(2,3,figsize=(12,9),sharex=True,sharey=True)

cax1 = fig.add_axes([0.039, 0.05, 0.305, 0.03]) 
cax2 = fig.add_axes([0.36, 0.05, 0.305, 0.03]) 
cax3 = fig.add_axes([0.68, 0.05, 0.305, 0.03]) 

a11 = ax11.contourf(param1['x_T']/1e3,param1['y_T']/1e3,h2mat(D1['hm'],param1),hlevs,cmap=cmocean.cm.balance)
fig.colorbar(a11,cax=cax1,orientation='horizontal',ticks=[-1.5,-1,-.5,0,.5,1,1.5])

a12 = ax12.contourf(param1['x_u']/1e3,param1['y_u']/1e3,u2mat(D1['um'],param1),ulevs,cmap=cmocean.cm.balance)
fig.colorbar(a12,cax=cax2,orientation='horizontal',ticks=[-.5,-.25,0,.25,.5])

a13 = ax13.contourf(param1['x_v']/1e3,param1['y_v']/1e3,v2mat(D1['vm'],param1),vlevs,cmap=cmocean.cm.balance)
fig.colorbar(a13,cax=cax3,orientation='horizontal',ticks=[-1,-.5,0,.5,1])

ax21.contourf(param2['x_T']/1e3,param2['y_T']/1e3,h2mat(D2['hm'],param2),hlevs,cmap=cmocean.cm.balance)        
ax22.contourf(param2['x_u']/1e3,param2['y_u']/1e3,u2mat(D2['um'],param2),ulevs,cmap=cmocean.cm.balance)    
ax23.contourf(param2['x_v']/1e3,param2['y_v']/1e3,v2mat(D2['vm'],param2),vlevs,cmap=cmocean.cm.balance)    

ax11.set_title(r'$\eta$ [m]')
ax12.set_title(r'$u$ [ms$^{-1}$]')
ax13.set_title(r'$v$ [ms$^{-1}$]')

ax11.set_yticks([])
ax21.set_yticks([])
ax21.set_xticks([])
ax22.set_xticks([])
ax23.set_xticks([])

ax11.set_ylabel(r'Low resolution, $\Delta x = 30$km')
ax21.set_ylabel(r'High resolution, $\Delta x = 7.5$km')

plt.tight_layout(rect=[0,.08,1,1])
plt.savefig(path+'compare/uvh_mean.png')
plt.close(fig)
