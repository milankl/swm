## PRODUCE MEAN PLOTS COMPARE FROM TWO RUNS
from __future__ import print_function

# path
import os
path = os.path.dirname(os.getcwd()) + '/'   # on level above
os.chdir(path)                              # change working directory

import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt
import cmocean

# OPTIONS
runfolder = [0,0]
print('Produce mean plots from run ' + str(runfolder))

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
etam_max = np.max(abs(D2['etam']))
etalevs = np.linspace(-etam_max,etam_max,64)

um_max = np.max(abs(D2['um']))
ulevs = np.linspace(-um_max,um_max,64)

vm_max = np.max(abs(D2['vm']))
vlevs = np.linspace(-vm_max,vm_max,64)

fig,axs = plt.subplots(2,3,figsize=(7,5),sharex=True,sharey=True)
plt.tight_layout(rect=[0,.06,1,0.99])
fig.subplots_adjust(wspace=0.03,hspace=0.03)
n = axs.shape[1]
caxs = [0,]*n
for i in range(n):
    pos = axs[-1,i].get_position()
    caxs[i] = fig.add_axes([pos.x0,0.095,pos.width,0.03])

a11 = axs[0,0].contourf(param1['x_T']/1e3,param1['y_T']/1e3,h2mat(D1['etam'],param1),etalevs,cmap=cmocean.cm.balance)
cb1 = fig.colorbar(a11,cax=caxs[0],orientation='horizontal',ticks=[-1.5,-1,-.5,0,.5,1,1.5])
cb1.set_label(r'[m]')

a12 = axs[0,1].contourf(param1['x_u']/1e3,param1['y_u']/1e3,u2mat(D1['um'],param1),ulevs,cmap=cmocean.cm.balance)
cb2 = fig.colorbar(a12,cax=caxs[1],orientation='horizontal',ticks=[-.5,-.25,0,.25,.5])
cb2.set_label(r'[ms$^{-1}$]')

a13 = axs[0,2].contourf(param1['x_v']/1e3,param1['y_v']/1e3,v2mat(D1['vm'],param1),vlevs,cmap=cmocean.cm.balance)
cb3 = fig.colorbar(a13,cax=caxs[2],orientation='horizontal',ticks=[-1,-.5,0,.5,1])   
cb3.set_label(r'[ms$^{-1}$]')

axs[1,0].contourf(param2['x_T']/1e3,param2['y_T']/1e3,h2mat(D2['etam'],param2),etalevs,cmap=cmocean.cm.balance)        
axs[1,1].contourf(param2['x_u']/1e3,param2['y_u']/1e3,u2mat(D2['um'],param2),ulevs,cmap=cmocean.cm.balance)    
axs[1,2].contourf(param2['x_v']/1e3,param2['y_v']/1e3,v2mat(D2['vm'],param2),vlevs,cmap=cmocean.cm.balance)    

axs[0,0].set_title(r'$\eta$')
axs[0,1].set_title(r'$u$')
axs[0,2].set_title(r'$v$')

axs[0,0].set_yticks([])
axs[0,0].set_xticks([])

axs[0,0].set_xlim(0,param1['Lx']/1e3)
axs[0,0].set_ylim(0,param1['Ly']/1e3)

axs[0,0].set_ylabel(r'LR, $\Delta x = 30$km')
axs[1,0].set_ylabel(r'HR, $\Delta x = 7.5$km')


plt.savefig(path+'figs/uvh_mean_2runs.png',dpi=150)
plt.close(fig)
