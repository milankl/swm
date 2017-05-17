## PLOT POWER MAPS FOR TWO RUNS
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
print('Plots for run ' + str(runfolder))

## read data

runpath1 = path+'data/run%04i' % runfolder[0]
D1 = np.load(runpath1+'/analysis/power_map.npy').all()
param1 = np.load(runpath1+'/param.npy').all()

runpath2 = path+'data/run%04i' % runfolder[1]
D2 = np.load(runpath2+'/analysis/power_map.npy').all()
param2 = np.load(runpath2+'/param.npy').all()

# functions
def h2mat(h,param):
    return h.reshape((param['ny'],param['nx']))

##
expo = 0.75   # stretching/squeezing for visualization

for v in ['InPower_T','ExPower_T']:
    D1[v] = np.sign(D1[v])*abs(D1[v])**expo
    D2[v] = np.sign(D2[v])*abs(D2[v])**expo

## PLOTTING   
fig,axs = plt.subplots(2,2,figsize=(5,5),sharex=True,sharey=True)
plt.tight_layout(rect=[0,.06,1,0.99])
fig.subplots_adjust(wspace=0.03,hspace=0.03)
pos = axs[-1,0].get_position()
pos2 = axs[-1,-1].get_position()
cax = fig.add_axes([pos.x0,0.09,pos2.x1-pos.x0,0.03])

levmax = abs(D2['InPower_T']).max()
levs = np.linspace(-levmax,levmax,64)
tiks = np.array([-0.06,-0.04,-0.02,0,0.02,0.04,0.06])

a11 = axs[0,0].contourf(param1['x_T']/1e3,param1['y_T']/1e3,h2mat(D1['InPower_T'],param1),levs,cmap=cm.balance,extend='both')
cb1 = fig.colorbar(a11,cax=cax,orientation='horizontal',ticks=np.sign(tiks)*abs(tiks)**expo)
cb1.set_ticklabels(tiks)
cb1.set_label(r'Power [Wm$^{-2}$]')

axs[0,1].contourf(param1['x_T']/1e3,param1['y_T']/1e3,h2mat(D1['ExPower_T'],param1),levs,cmap=cm.balance,extend='both')

axs[1,0].contourf(param2['x_T']/1e3,param2['y_T']/1e3,h2mat(D2['InPower_T'],param2),levs,cmap=cm.balance,extend='both')        
axs[1,1].contourf(param2['x_T']/1e3,param2['y_T']/1e3,h2mat(D2['ExPower_T'],param2),levs,cmap=cm.balance,extend='both')      

axs[0,0].set_title('wind forcing')
axs[0,1].set_title('lateral mixing')

axs[0,0].set_yticks([])
axs[0,0].set_xticks([])

axs[0,0].set_xlim(0,param1['Lx']/1e3)
axs[0,0].set_ylim(0,param1['Ly']/1e3)

axs[0,0].set_ylabel(r'LR, $\Delta x = 30$km')
axs[1,0].set_ylabel(r'HR, $\Delta x = 7.5$km')

plt.savefig(path+'figs/power_maps.png',dpi=150)
plt.close(fig)