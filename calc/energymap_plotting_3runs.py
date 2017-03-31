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
from cmocean import cm

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

runpath3 = path+'stoch/data/run%04i' % 12
D3 = np.load(runpath3+'/analysis/mean.npy').all()
param3 = np.load(runpath3+'/param.npy').all()

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

expo = 1.

for v in ['mke','eke']:
    D1[v] = (D1[v]/param1['rho']/param1['H'])**expo
    D2[v] = (D2[v]/param1['rho']/param1['H'])**expo
    D3[v] = (D3[v]/param1['rho']/param1['H'])**expo

for v in ['epe','mpe']:
    D1[v] = (D1[v]/param1['rho']/param1['g'])**expo
    D2[v] = (D2[v]/param1['rho']/param1['g'])**expo
    D3[v] = (D3[v]/param1['rho']/param1['g'])**expo

## PLOTTING   
levs=[]
for v in ['mke','eke','mpe','epe']:
    maxs = np.max((np.max(D1[v]),np.max(D2[v])))
    #mins = np.min((np.min(D1[v]),np.min(D2[v])))
    levs.append(np.linspace(0,maxs*0.95,64))

fig,axs = plt.subplots(3,4,figsize=(12,9),sharex=True,sharey=True)
plt.tight_layout(rect=[0,.08,1,0.98])
fig.subplots_adjust(wspace=0.03,hspace=0.03)
n = axs.shape[1]
caxs = [0,]*n
for i in range(n):
    pos = axs[-1,i].get_position()
    caxs[i] = fig.add_axes([pos.x0,0.07,pos.width,0.03])

qaxs = np.empty_like(axs)

tiks = [np.array([0,.25,.5,.75]),np.array([0,1,2,3]),np.array([0,0.25,.5,0.75,1]),np.array([0,1,2,3,4,5])]

qaxs[0,0] = axs[0,0].contourf(param1['x_T'],param1['y_T'],h2mat(D1['mke'],param1),levs[0],cmap=cm.thermal,extend='max')
cb1 = fig.colorbar(qaxs[0,0],cax=caxs[0],orientation='horizontal',ticks=tiks[0]**expo)
cb1.set_ticklabels(tiks[0])
cb1.set_label(r'[m$^2$s$^{-2}$]')

qaxs[0,1] = axs[0,1].contourf(param1['x_T'],param1['y_T'],h2mat(D1['eke'],param1),levs[1],cmap=cm.thermal,extend='max')
cb2 = fig.colorbar(qaxs[0,1],cax=caxs[1],orientation='horizontal',ticks=tiks[1]**expo)
cb2.set_ticklabels(tiks[1].astype(np.int))
cb2.set_label(r'[m$^2$s$^{-2}$]')

qaxs[0,2] = axs[0,2].contourf(param1['x_T'],param1['y_T'],h2mat(D1['mpe'],param1),levs[2],cmap=cm.thermal,extend='max')
cb3 = fig.colorbar(qaxs[0,2],cax=caxs[2],orientation='horizontal',ticks=tiks[2]**expo)
cb3.set_ticklabels(tiks[2])
cb3.set_label(r'[m$^2$]')

qaxs[0,3] = axs[0,3].contourf(param1['x_T'],param1['y_T'],h2mat(D1['epe'],param1),levs[3],cmap=cm.thermal,extend='max')
cb4 = fig.colorbar(qaxs[0,3],cax=caxs[3],orientation='horizontal',ticks=tiks[3]**expo)
cb4.set_ticklabels(tiks[3])
cb4.set_label(r'[m$^2$]')

axs[1,0].contourf(param3['x_T'],param3['y_T'],h2mat(D3['mke'],param3),levs[0],cmap=cm.thermal,extend='max')
axs[1,1].contourf(param3['x_T'],param3['y_T'],h2mat(D3['eke'],param3),levs[1],cmap=cm.thermal,extend='max')    
axs[1,2].contourf(param3['x_T'],param3['y_T'],h2mat(D3['mpe'],param3),levs[2],cmap=cm.thermal,extend='max')    
axs[1,3].contourf(param3['x_T'],param3['y_T'],h2mat(D3['epe'],param3),levs[3],cmap=cm.thermal,extend='max')    

axs[2,0].contourf(param2['x_T'],param2['y_T'],h2mat(D2['mke'],param2),levs[0],cmap=cm.thermal,extend='max')
axs[2,1].contourf(param2['x_T'],param2['y_T'],h2mat(D2['eke'],param2),levs[1],cmap=cm.thermal,extend='max')    
axs[2,2].contourf(param2['x_T'],param2['y_T'],h2mat(D2['mpe'],param2),levs[2],cmap=cm.thermal,extend='max')    
axs[2,3].contourf(param2['x_T'],param2['y_T'],h2mat(D2['epe'],param2),levs[3],cmap=cm.thermal,extend='max')    

axs[0,0].set_title('MKE')
axs[0,1].set_title('EKE')
axs[0,2].set_title('MPE')
axs[0,3].set_title('EPE')

axs[0,0].set_xticks([])
axs[0,0].set_yticks([])

axs[0,0].set_xlim(0,param1['Lx'])
axs[0,0].set_ylim(0,param1['Ly'])

axs[0,0].set_ylabel(r'Low resolution, $\Delta x = 30$km')
axs[1,0].set_ylabel(r'Low resolution + backscatter')
axs[2,0].set_ylabel(r'High resolution, $\Delta x = 7.5$km')


plt.savefig(path+'compare/Energy_mean.png')
#plt.close(fig)
#plt.show()
