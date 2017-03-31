## RUN SHALLOW WATER MODEL
from __future__ import print_function
path = '/home/mkloewer/python/swm/'

# import modules
import os; os.chdir(path)
import numpy as np
from scipy import sparse
from scipy.interpolate import RegularGridInterpolator as RIG
import time as tictoc
from netCDF4 import Dataset
import glob
import zipfile

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as manimation

FFMpegWriter = manimation.writers['ffmpeg']
writer = FFMpegWriter(fps=25)

## import all functions
exec(open(path+'swm_param.py').read())
exec(open(path+'swm_operators.py').read())
exec(open(path+'swm_rhs.py').read())
exec(open(path+'swm_integration.py').read())
exec(open(path+'swm_output.py').read())

## read image
from matplotlib.image import imread
img = imread('analyses/floats/richyrich.png').swapaxes(0,1)[:,:,:]
nx,ny,_ = img.shape
img = img.reshape((-1,3))

## read data
runfolder = 10
runpath = path+'data/run%04i' % runfolder

tmax = 500

ncu = Dataset(runpath+'/u.nc')
u = ncu['u'][:500,:,:]
ncu.close()
print('u read.')

ncv = Dataset(runpath+'/v.nc')
v = ncv['v'][:500,:,:]
ncv.close()
print('v read.')

print('netCDF data read.')

# read param
global param
param = np.load(runpath+'/param.npy').all()
param['dat_type'] = np.float32

# import functions
param['output'] = 0

set_grad_mat()
set_interp_mat()
set_lapl_mat()
set_coriolis()

tlen = len(time) 

ntime = u.shape[0]
dt = time[1]-time[0]

## final conditions (='initial conditions')
xf,yf = np.meshgrid(np.linspace(param['x_T'][150],param['x_T'][350],nx),np.linspace(param['y_T'][350],param['y_T'][150],ny))
xf = xf.T.flatten()
yf = yf.T.flatten()

# renaming for convenience
x_u = param['x_u']
y_u = param['y_u']

x_v = param['x_v']
y_v = param['y_v']

# preallocate
X = np.empty((ntime,xf.shape[0]))
Y = np.empty((ntime,xf.shape[0]))

# final conditions
X[-1,:] = xf
Y[-1,:] = yf

for i in range(ntime-1,0,-1):
    Iu = RIG((x_u,y_u),u[i,:,:].T,bounds_error=False,fill_value=0)
    Iv = RIG((x_v,y_v),v[i,:,:].T,bounds_error=False,fill_value=0)

    X[i-1,:] = X[i,:] - dt*Iu((X[i,:],Y[i,:]))
    Y[i-1,:] = Y[i,:] - dt*Iv((X[i,:],Y[i,:]))
    if i % 10 == 0:
        print(i)
    
##
fig,ax = plt.subplots(1,1)

ax.set_xlim(0,param['Lx'])
ax.set_ylim(0,param['Lx'])
ax.set_xticks([])
ax.set_yticks([])
plt.tight_layout()
scat = ax.scatter(X[0,:],Y[0,:],10,facecolors=img,edgecolors='none',marker='s')

def update(j):
    scat.set_offsets(np.vstack((X[j,:],Y[j,:])).T)
    
with writer.saving(fig, "analyses/floats/richyrich.avi", 200):
    for j in list(range(tmax))+[tmax-1]*50:
        update(j)
        writer.grab_frame()