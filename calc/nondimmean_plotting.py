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
##
for v in ['Re','Ro','Ek']:
    D1[v] = np.log10(D1[v])
    D2[v] = np.log10(D2[v])

## PLOTTING   
levs = [0,]*3
#for v in ['Re','Ro','Ek']:
#    maxs = np.max((np.max(D1[v]),np.max(D2[v])))
#    mins = np.min((np.min(D1[v]),np.min(D2[v])))
#    levs.append(np.linspace(mins,maxs,9))
levs[0] = np.arange(0,3.25,.25)
levs[1] = np.arange(-3,0,.25)
levs[2] = np.arange(-5,-2,.25)


fig,((ax11,ax12,ax13),(ax21,ax22,ax23)) = plt.subplots(2,3,figsize=(12,9),sharex=True,sharey=True)

cax1 = fig.add_axes([0.039, 0.05, 0.305, 0.03]) 
cax2 = fig.add_axes([0.36, 0.05, 0.305, 0.03]) 
cax3 = fig.add_axes([0.68, 0.05, 0.305, 0.03]) 

a11 = ax11.contourf(param1['x_T']/1e3,param1['y_T']/1e3,h2mat(D1['Re'],param1),levs[0],cmap='viridis',extend='min')
fig.colorbar(a11,cax=cax1,orientation='horizontal',ticks=[0,.5,1,1.5,2,2.5])

a12 = ax12.contourf(param1['x_T']/1e3,param1['y_T']/1e3,h2mat(D1['Ro'],param1),levs[1],cmap='viridis',extend='min')
fig.colorbar(a12,cax=cax2,orientation='horizontal',ticks=[-3,-2.5,-2,-1.5,-1,-0.5])

a13 = ax13.contourf(param1['x_T']/1e3,param1['y_T']/1e3,h2mat(D1['Ek'],param1),levs[2],cmap='viridis',extend='max')
fig.colorbar(a13,cax=cax3,orientation='horizontal',ticks=[-5,-4.5,-4,-3.5,-3,-2.5])

ax21.contourf(param2['x_T']/1e3,param2['y_T']/1e3,h2mat(D2['Re'],param2),levs[0],cmap='viridis',extend='min')
ax22.contourf(param2['x_T']/1e3,param2['y_T']/1e3,h2mat(D2['Ro'],param2),levs[1],cmap='viridis',extend='min')    
ax23.contourf(param2['x_T']/1e3,param2['y_T']/1e3,h2mat(D2['Ek'],param2),levs[2],cmap='viridis',extend='max')    

ax11.set_title(r'log$_{10}(Re)$')
ax12.set_title(r'log$_{10}(Ro)$')
ax13.set_title(r'log$_{10}(Ek)$')

ax11.set_yticks([])
ax21.set_yticks([])
ax21.set_xticks([])
ax22.set_xticks([])
ax23.set_xticks([])

ax11.set_ylabel(r'Low resolution, $\Delta x = 30$km')
ax21.set_ylabel(r'High resolution, $\Delta x = 7.5$km')

plt.tight_layout(rect=[0,.08,1,1])
plt.savefig(path+'compare/ReRoEk_mean.png')
plt.close(fig)
