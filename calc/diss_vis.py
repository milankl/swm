## UNDERSTANDING DISSIPATION
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
from cmocean import cm

# import functions
exec(open(path+'swm_param.py').read())
exec(open(path+'swm_operators.py').read())
exec(open(path+'swm_output.py').read())
exec(open(path+'stoch/swm_rhs.py').read())

## OPTIONS
runfolder = [3,10]

## read data
runpath = path+'data/run%04i' % runfolder[0]
ncu = Dataset(runpath+'/u.nc')
ncv = Dataset(runpath+'/v.nc')
nch = Dataset(runpath+'/h.nc')

u1 = ncu['u'][-1000,:,:].flatten()
v1 = ncv['v'][-1000,:,:].flatten()
h1 = nch['h'][-1000,:,:].flatten()
print('netCDF data read.')

# close netcdfs
ncu.close()
ncv.close()
nch.close()

param1 = np.load(runpath+'/param.npy').all()
param1['output'] = 0

##
runpath = path+'data/run%04i' % runfolder[1]
ncu = Dataset(runpath+'/u.nc')
ncv = Dataset(runpath+'/v.nc')
nch = Dataset(runpath+'/h.nc')

u2 = ncu['u'][-1,:,:].flatten()
v2 = ncv['v'][-1,:,:].flatten()
h2 = nch['h'][-1,:,:].flatten()
print('netCDF data read.')

# close netcdfs
ncu.close()
ncv.close()
nch.close()

param2 = np.load(runpath+'/param.npy').all()
param2['output'] = 0

##
global param

def Q1(u,v,h):

    set_grad_mat()
    set_interp_mat()
    set_lapl_mat()
    set_coriolis()
    
    h_u = ITu.dot(h)    # h on u-grid
    h_v = ITv.dot(h)    # h on v-grid    
    h_q = ITq.dot(h)    # h on q-grid
    
    dudx = Gux.dot(u)
    dudy = Guy.dot(u)
    dvdx = Gvx.dot(v)
    dvdy = Gvy.dot(v)
    
    #diff_u, diff_v = mixing(u,v,h,h_q,h_u,h_v)              # harmonic
    #bidiff_u, bidiff_v = mixing(diff_u,diff_v,h,h_q,h_u,h_v)    # apply twice = biharmonic
    
    diss_e = (dudx*Gux.dot(Lu.dot(u)) + dvdy*Gvy.dot(Lv.dot(v)))
    
    return diss_e
    
def Q2(u,v,h):

    set_grad_mat()
    set_interp_mat()
    set_lapl_mat()
    set_coriolis()
    
    h_u = ITu.dot(h)    # h on u-grid
    h_v = ITv.dot(h)    # h on v-grid    
    h_q = ITq.dot(h)    # h on q-grid
    
    dudx = Gux.dot(u)
    dudy = Guy.dot(u)
    dvdx = Gvx.dot(v)
    dvdy = Gvy.dot(v)
    
    diff_u, diff_v = mixing(u,v,h,h_q,h_u,h_v)              # harmonic
    bidiff_u, bidiff_v = mixing(diff_u,diff_v,h,h_q,h_u,h_v)    # apply twice = biharmonic
    
    #bidiff_u = LLu.dot(u)
    #bidiff_v = LLv.dot(v)
    
    #diss_e = param['B']*(dudx*Gux.dot(Lu.dot(u)) + dvdy*Gvy.dot(Lv.dot(v)))
    
    return IuT.dot(u*bidiff_u) + IvT.dot(v*bidiff_v)

param = param2
q = Q1(u2,v2,h2+param['H'])
q2 = Q2(u2,v2,h2+param['H'])

## plotting
fig,(ax1,ax2) = plt.subplots(1,2,sharex=True,sharey=True,figsize=(12,6))
plt.tight_layout(rect=[0,0,1,0.95])

pq = q*2e16*100
pq2 = q2*2e16

levs = np.linspace(-abs(pq2).max()*0.02,abs(pq2).max()*0.02,64)

param = param2
c = ax1.contourf(param2['x_T'],param2['y_T'],h2mat(pq),levs,extend='both',cmap='RdBu_r')
c = ax2.contourf(param2['x_T'],param2['y_T'],h2mat(pq2),levs,extend='both',cmap='RdBu_r')

ax1.set_xticks([])
ax1.set_yticks([])
ax1.set_xlabel('x')
ax2.set_xlabel('x')
ax1.set_ylabel('y')

ax1.set_xlim(571e3,1824e3)
ax1.set_ylim(1051e3,2131e3)

ax1.set_title(r'$100*\nabla\mathbf{u} \cdot \nabla(\nabla^2\mathbf{u})$')
ax2.set_title(r'$\mathbf{u} \cdot \nabla^4\mathbf{u}$')

cbar = plt.colorbar(c,ax=(ax1,ax2),ticks=[-1,-0.5,0,0.5,1])
cbar.set_label(r'[$2 \cdot 10^{-16}m^{-2}s^{-2}$]')
plt.show()