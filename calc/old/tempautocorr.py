## COMPUTE AND PRODUCE TIMESCALE PLOTS
from __future__ import print_function
path = '/home/mkloewer/python/swm/'
import os; os.chdir(path) # change working directory
import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt
import time as tictoc
from netCDF4 import Dataset
import glob
import cmocean
from matplotlib.colors import BoundaryNorm

## OPTIONS
runfolder = 2
print(('Running temporal autocorrelation scales from run %i') % runfolder)

## read data
runpath = path+'data/nesh/run%04i' % runfolder
ncu = Dataset(runpath+'/u.nc')
ncv = Dataset(runpath+'/v.nc')
nch = Dataset(runpath+'/h.nc')

u = ncu['u'][:,:,:]
v = ncv['v'][:,:,:]
h = nch['h'][:,:,:]
time = nch['t'][:]   # in seconds
t = time / 3600. / 24.  # in days
print('netCDF data read.')

# close netcdfs
ncu.close()
ncv.close()
nch.close()

# read param
global param
param = np.load(runpath+'/param.npy').all()
tlen = time.shape[0]

# import functions
exec(open(path+'swm_param.py').read())
exec(open(path+'swm_operators.py').read())
exec(open(path+'swm_output.py').read())
param['output'] = 0

set_grad_mat()
set_interp_mat()
set_lapl_mat()
set_coriolis()

tlen = len(time)
##
u = u.reshape((tlen,param['Nu'])).T
v = v.reshape((tlen,param['Nv'])).T
h = h.reshape((tlen,param['NT'])).T
z = Gvx.dot(v) - Guy.dot(u)

def acfast_parallel(A,lagmax):
    """ Fast parallel autocorrelation, A is supposed to be of shape time x variable. Mean and standard deviation of variables in A are assumed to be constant over time. Only appropriate for lagmax << len(time)."""
    n,p = A.shape
    A = A - A.mean(axis=0)   # anomalies
    # A but normalised and standardised on the time axis
    A = A / np.sqrt((A**2).sum(axis=0) / n)  
    tau = np.empty((p,lagmax))  # preallocate
    tau[:,0] = 1    # do not compute for lag 0
    for l in range(1,lagmax):
        tau[:,l] = (A[l:,:]*A[:-l,:]).mean(axis=0)

    return tau

def find_crossing_parallel(A,a,dt):
    """ A = var x time. crossing with a, only first occurence. """
    p,n = A.shape
    all_ind_var,all_ind_time = np.where(np.diff(np.signbit(A-a))) # all crossings as indices
    
    # preallocate time array when crossings occur, zero hence meaning no crossing
    t = np.zeros(p,dtype=int)    
    # trick: start from the last ([::-1]-operation, to neglect all but the first occurence
    t[all_ind_var[::-1]] = all_ind_time[::-1]  # returning this is already an estimate for the crossing
    
    # linear interpolation
    v = np.arange(p)
    A1,A0 = A[v,t+1],A[v,t]
    return dt*(t+1-(a-A1)/(A0 - A1))

def decorr_tscale(quant):
    a = 1/np.exp(1) # threshold
    lagmax = 100
    dt = t[1] - t[0]
    acf = acfast_parallel(quant.T,lagmax)
    tau = find_crossing_parallel(acf,a,dt)
    return np.ma.masked_array(tau,mask=np.isnan(tau))
    
htau = decorr_tscale(h)
print('Timescales of h done.')
utau = decorr_tscale(u)
print('Timescales of u done.')
vtau = decorr_tscale(v)
print('Timescales of v done.')
ztau = decorr_tscale(z)
print('Timescales of z done.')
    
## plot
p = 1

# H
fig,ax = plt.subplots(1,1,figsize=(12,9))

levs = np.linspace(np.percentile(htau.compressed(),p),np.percentile(htau.compressed(),100-p),64)

c = ax.contourf(param['x_T']/1e3,param['y_T']/1e3,h2mat(htau),levs,extend='both',cmap='viridis') 
cb = plt.colorbar(c,ax=ax)
cb.set_label(r'$\tau$ [days]')

ax.set_ylabel('y [km]')
ax.set_xlabel('x [km]')
ax.set_title('Decorrelation time scale: $h$')

plt.tight_layout()
plt.savefig(runpath+'/plots/h_tscale.png')
plt.close(fig)

# U VELOCITY
fig,ax = plt.subplots(1,1,figsize=(12,9))

levs = np.linspace(np.percentile(utau.compressed(),p),np.percentile(utau.compressed(),100-p),64)

c = ax.contourf(param['x_u']/1e3,param['y_u']/1e3,u2mat(utau),levs,extend='both',cmap='viridis') 
cb = plt.colorbar(c,ax=ax)
cb.set_label(r'$\tau$ [days]')

ax.set_ylabel('y [km]')
ax.set_xlabel('x [km]')
ax.set_title('Decorrelation time scale: $u$')

plt.tight_layout()
plt.savefig(runpath+'/plots/u_tscale.png')
plt.close(fig)

# V VELOCITY
fig,ax = plt.subplots(1,1,figsize=(12,9))

levs = np.linspace(np.percentile(vtau.compressed(),p),np.percentile(vtau.compressed(),100-p),64)

c = ax.contourf(param['x_v']/1e3,param['y_v']/1e3,v2mat(vtau),levs,extend='both',cmap='viridis') 
cb = plt.colorbar(c,ax=ax)
cb.set_label(r'$\tau$ [days]')

ax.set_ylabel('y [km]')
ax.set_xlabel('x [km]')
ax.set_title('Decorrelation time scale: $v$')

plt.tight_layout()
plt.savefig(runpath+'/plots/v_tscale.png')
plt.close(fig)

# RELATIVE VORTICITY
fig,ax = plt.subplots(1,1,figsize=(12,9))

levs = np.linspace(np.percentile(ztau.compressed(),p),np.percentile(ztau.compressed(),100-p),64)

c = ax.contourf(param['x_q']/1e3,param['y_q']/1e3,q2mat(ztau),levs,extend='both',cmap='viridis') 
cb = plt.colorbar(c,ax=ax)
cb.set_label(r'$\tau$ [days]')

ax.set_ylabel('y [km]')
ax.set_xlabel('x [km]')
ax.set_title('Decorrelation time scale: Relative vorticity')

plt.tight_layout()
plt.savefig(runpath+'/plots/rel_vort_tscale.png')
plt.close(fig)
