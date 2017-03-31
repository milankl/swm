## COMPUTE AND PRODUCE TWO SUBPLOTS WITH RELAITVE VORTICITY FROM DIFFERENT RUNS
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
import cmocean
from matplotlib.colors import BoundaryNorm

## OPTIONS
runfolder = [20,21,19]

## read data
runpath = path+'data/run%04i' % runfolder[0]
ncu = Dataset(runpath+'/u.nc')
ncv = Dataset(runpath+'/v.nc')

t = 40

u1 = ncu['u'][t,:,:].flatten()
v1 = ncv['v'][t,:,:].flatten()
print('netCDF data read.')

# close netcdfs
ncu.close()
ncv.close()

# read param
global param
param = np.load(runpath+'/param.npy').all()

# import functions
exec(open(path+'swm_param.py').read())
exec(open(path+'swm_operators.py').read())
exec(open(path+'swm_output.py').read())
param['output'] = 0

set_grad_mat()
set_interp_mat()
set_lapl_mat()
set_coriolis()

z1 = Gvx.dot(v1) - Guy.dot(u1)

## SECOND FILE
runpath = path+'data/run%04i' % runfolder[1]
ncu = Dataset(runpath+'/u.nc')
ncv = Dataset(runpath+'/v.nc')

u2 = ncu['u'][t,:,:].flatten()
v2 = ncv['v'][t,:,:].flatten()
print('netCDF data read.')

# close netcdfs
ncu.close()
ncv.close()

# read param
param = np.load(runpath+'/param.npy').all()

# import functions
exec(open(path+'swm_param.py').read())
exec(open(path+'swm_operators.py').read())
exec(open(path+'swm_output.py').read())
param['output'] = 0

set_grad_mat()
set_interp_mat()
set_lapl_mat()
set_coriolis()

z2 = Gvx.dot(v2) - Guy.dot(u2)

## THIRD FILE
runpath = path+'data/run%04i' % runfolder[2]
ncu = Dataset(runpath+'/u.nc')
ncv = Dataset(runpath+'/v.nc')

u3 = ncu['u'][t,:,:].flatten()
v3 = ncv['v'][t,:,:].flatten()
print('netCDF data read.')

# close netcdfs
ncu.close()
ncv.close()

# read param
param = np.load(runpath+'/param.npy').all()

# import functions
exec(open(path+'swm_param.py').read())
exec(open(path+'swm_operators.py').read())
exec(open(path+'swm_output.py').read())
param['output'] = 0

set_grad_mat()
set_interp_mat()
set_lapl_mat()
set_coriolis()

z3 = Gvx.dot(v3) - Guy.dot(u3)

## REL VORTICITY FINAL
z1 = np.sign(z1)*abs(z1)**(1/2.)
z2 = np.sign(z2)*abs(z2)**(1/2.)
z3 = np.sign(z3)*abs(z3)**(1/2.)

z1_max = np.percentile(abs(z1),98)
levs = np.linspace(-z1_max,z1_max,64)

##
fig,(ax1,ax2) = plt.subplots(1,2,figsize=(14,6),sharex=True,sharey=True)
ax1.contourf(param['x_q']/1e3,param['y_q']/1e3,q2mat(z1),levs,cmap='RdBu_r',extend='both')    

ax1.set_xlim(0,param['Lx']/1e3)
ax1.set_ylim(0,param['Ly']/1e3)
ax1.set_xticks([])
ax1.set_yticks([])
#ax1.set_title(r'(a) Relative vorticity $\zeta_a$, precision: 64bit')
ax1.text(0.96,0.96,r'(a)',ha='right',va='top',fontsize=20,transform=ax1.transAxes,bbox=dict(fc='w',ec='k'))

ax2.contourf(param['x_q']/1e3,param['y_q']/1e3,q2mat(z2),levs,cmap='RdBu_r',extend='both')
#ax2.set_title(r'(b) Relative vorticity $\zeta_b$, precision: 32bit')
ax2.text(0.96,0.96,r'(b)',ha='right',va='top',fontsize=20,transform=ax2.transAxes,bbox=dict(fc='w',ec='k'))

plt.tight_layout()
#plt.savefig(runpath+'/plots/relvort_final.png')
#plt.close(fig)
plt.show()