from __future__ import print_function

path = '/network/home/aopp/kloewer/strix/'
dpath = '/network/aopp/cirrus/pred/kloewer/swm_back_ronew/'
outpath = '/network/home/aopp/kloewer/swm/paperplot/'

import os; os.chdir(path) # change working directory
import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt
import time as tictoc
from netCDF4 import Dataset
import glob
from matplotlib.colors import BoundaryNorm,LogNorm
import cmocean

plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'

# OPTIONS
runfolder = [0,6,2]
print('Plots for run ' + str(runfolder))

## read data

runpath1 = path+'run%04i' % runfolder[0]
D1 = np.load(runpath1+'/analysis/power_map.npy').all()
param1 = np.load(runpath1+'/param.npy').all()

runpath2 = path+'run%04i' % runfolder[1]
D2 = np.load(runpath2+'/analysis/power_map.npy').all()
param2 = np.load(runpath2+'/param.npy').all()

runpath3 = dpath+'run%04i' % runfolder[2]
D3 = np.load(runpath3+'/analysis/power_map.npy').all()
param3 = np.load(runpath3+'/param.npy').all()
print(param3['n_diss'])

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

m = [0]*(len(runfolder)*3) 
s = 1e3

# LOW RESOLUTION RUN
in1 = h2mat(D1['InPower_T'],param1)*s
m[0] = in1.mean()
in1 = np.sign(in1)*np.sqrt(abs(in1))

bf1 = h2mat(D1['BfricPower_T'],param1)*s
m[1] = bf1.mean()
bf1 = np.sign(bf1)*np.sqrt(abs(bf1))

ex1 = h2mat(D1['ExPower_T'],param1)*s
m[2] = ex1.mean()
ex1 = np.sign(ex1)*np.sqrt(abs(ex1))

# BACKSCATTER RUN
in3 = h2mat(D3['InPower_T'],param3)*s
m[3] = in3.mean()
in3 = np.sign(in3)*np.sqrt(abs(in3))

bf3 = h2mat(D3['BfricPower_T'],param3)*s
m[4] = bf3.mean()
bf3 = np.sign(bf3)*np.sqrt(abs(bf3))

D3['ExPower_T'] += D3['BackPower_T']

ex3 = h2mat(D3['ExPower_T'],param3)*s
m[5] = ex3.mean()
ex3 = np.sign(ex3)*np.sqrt(abs(ex3))

# bs = h2mat(D3['BackPower_T'],param3)*s
# m[6] = bs.mean()
# bs = np.sign(bs)*np.sqrt(abs(bs))

# HIGH RESOLUTION RUN
in2 = h2mat(D2['InPower_T'],param2)*s
m[6] = in2.mean()
in2 = np.sign(in2)*np.sqrt(abs(in2))

bf2 = h2mat(D2['BfricPower_T'],param2)*s
m[7] = bf2.mean()
bf2 = np.sign(bf2)*np.sqrt(abs(bf2))

ex2 = h2mat(D2['ExPower_T'],param2)*s
m[8] = ex2.mean()
ex2 = np.sign(ex2)*np.sqrt(abs(ex2))

mround = [np.round(mi,2) for mi in m]

## PLOTTING   

fig,axs = plt.subplots(3,3,figsize=(10,9),sharex=True,sharey=True)

plt.tight_layout(rect=[0,.07,1,0.97])
fig.subplots_adjust(wspace=0.03,hspace=0.03)

pos = axs[-1,0].get_position()
pos2 = axs[-1,-1].get_position()
cax = fig.add_axes([pos.x0,0.06,pos2.x1-pos.x0,0.02])

levs = np.linspace(-np.sqrt(90),np.sqrt(90),64)
tik = np.array([-90,-60,-30,-10,-3,0,3,10,30,60,90])
tik = np.sign(tik)*np.sqrt(abs(tik))

q1 = axs[0,0].contourf(param1['x_T'],param1['y_T'],in1,levs,cmap=cmocean.cm.balance,extend='both')
cbar = fig.colorbar(q1,cax=cax,orientation='horizontal',ticks=tik)
cbar.set_label(r'Power [Wm$^{-2} \cdot 10^{-3}$]')
cbar.set_ticklabels(np.round(tik**2*np.sign(tik)).astype(int))

axs[1,0].contourf(param3['x_T'],param3['y_T'],in3,levs,cmap=cmocean.cm.balance,extend='both')    
axs[2,0].contourf(param2['x_T'],param2['y_T'],in2,levs,cmap=cmocean.cm.balance,extend='both')    
    
axs[0,2].contourf(param1['x_T'],param1['y_T'],ex1,levs,cmap=cmocean.cm.balance,extend='both')
axs[1,2].contourf(param3['x_T'],param3['y_T'],ex3,levs,cmap=cmocean.cm.balance,extend='both')
axs[2,2].contourf(param2['x_T'],param2['y_T'],ex2,levs,cmap=cmocean.cm.balance,extend='both')  

axs[0,1].contourf(param1['x_T'],param1['y_T'],bf1,levs,cmap=cmocean.cm.balance,extend='both')
axs[1,1].contourf(param3['x_T'],param3['y_T'],bf3,levs,cmap=cmocean.cm.balance,extend='both')
axs[2,1].contourf(param2['x_T'],param2['y_T'],bf2,levs,cmap=cmocean.cm.balance,extend='both')

#axs[1,3].contourf(param3['x_T'],param3['y_T'],bs,levs,cmap=cmocean.cm.balance,extend='both')


axs[0,0].set_title('Wind forcing power',loc='left')
axs[0,1].set_title('Bottom friction power',loc='left')
axs[0,2].set_title('Biharmonic viscosity\n+ backscatter power',loc='left')
#axs[1,3].set_title('Backscatter power',loc='left')

abc = 'abcdefghijkl'
abci = 0
for i,axcol in enumerate(axs):
    for j,ax in enumerate(axcol):
#        if not(i in [0,2] and j == 3):
        plt.text(0.93,0.93,abc[abci],transform=ax.transAxes,fontweight='bold')
        plt.text(0.965,0.87,"%.2f" % mround[abci],transform=ax.transAxes,ha='right')
        abci += 1

axs[0,0].set_xticks([])
axs[0,0].set_yticks([])
axs[0,0].set_xlim(0,param1['Lx'])
axs[0,0].set_ylim(0,param1['Ly'])

axs[0,0].set_ylabel(r'Low resolution, $\Delta x = $30km')
axs[1,0].set_ylabel(r'LR + moderate backscatter')
axs[2,0].set_ylabel(r'High resolution, $\Delta x = $7.5km')

axs[-1,0].set_xlabel(r'$x$')
axs[-1,1].set_xlabel(r'$x$')
axs[-1,2].set_xlabel(r'$x$')
#axs[-1,3].set_xlabel(r'$x$')

#axs[0,3].set_frame_on(False)
#axs[2,3].set_frame_on(False)

axs[0,2].set_ylabel(r'$y$')
axs[0,2].yaxis.set_label_position('right')

axs[1,2].set_ylabel(r'$y$')
axs[1,2].yaxis.set_label_position('right')

axs[2,2].set_ylabel(r'$y$')
axs[2,2].yaxis.set_label_position('right')

plt.savefig(outpath+'plots/power_maps_bs.png',dpi=150)
plt.close(fig)
