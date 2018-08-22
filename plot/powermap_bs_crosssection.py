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
#ex1 = np.sign(ex1)*np.sqrt(abs(ex1))

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
#ex3 = np.sign(ex3)*np.sqrt(abs(ex3))

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
#ex2 = np.sign(ex2)*np.sqrt(abs(ex2))

mround = [np.round(mi,2) for mi in m]

## PLOTTING   

fig,(ax,ax2) = plt.subplots(2,1)

s1 = slice(32,96,1)     # for low resolution from y/4 to 3y/4
s2 = slice(128,384,1)   # for high resolution

ax.plot((param1['x_T']-param1['dx']/2.)/1e3,ex1[s1,:].mean(axis=0),'C0',drawstyle='steps-post')
ax.plot((param2['x_T']-param2['dx']/2.)/1e3,ex2[s2,:].mean(axis=0),'C2',drawstyle='steps-post')
ax.plot((param3['x_T']-param3['dx']/2.)/1e3,ex3[s1,:].mean(axis=0),'C3',drawstyle='steps-post')

ax2.plot((param1['x_T']-param1['dx']/2.)/1e3,ex1[s1,:].mean(axis=0),'C0',label='Low resolution, $\Delta x = 30$km',drawstyle='steps-post')
ax2.plot((param2['x_T']-param2['dx']/2.)/1e3,ex2[s2,:].mean(axis=0),'C2',label='High resolution, $\Delta x = 7.5$km',drawstyle='steps-post')
ax2.plot((param3['x_T']-param3['dx']/2.)/1e3,ex3[s1,:].mean(axis=0),'C3',label='LR + moderate backscatter',drawstyle='steps-post')

ax.set_ylim(-20,20)
ax.set_xlim(0,param1['Lx']/1e3)

#ax2.set_ylim(-20,20)
ax2.set_xlim(0,param1['Lx']/1e3/6)


ax.set_title(r'Biharm. viscosity + backscatter power, $y=\frac{L_y}{4} \colon \frac{3L_y}{4}$@average')
ax2.set_title('Zoom in in $x$, zoom out in y-axis')
ax.set_title('(a)',loc='left')
ax2.set_title('(b)',loc='left')


ax.set_ylabel(r'mW/m$^2$')
ax2.set_ylabel(r'mW/m$^2$')
ax.set_xlabel('$x$ [km]')
ax2.set_xlabel('$x$ [km]')

ax2.legend(loc=4)

plt.tight_layout()
plt.savefig(outpath+'plots/power_maps_bs_crosssection.png',dpi=150)
plt.close(fig)
