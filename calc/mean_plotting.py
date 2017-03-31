## PRODUCE MEAN PLOTS and compute
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

# OPTIONS
runfolder = 13
print(('Creating mean plots from run %i') % runfolder)

## read data
runpath = path+'data/run%04i' % runfolder

D = np.load(runpath+'/analysis/mean.npy').all()
for k in list(D.keys()):
   exec(k+' = D["'+k+'"]')
   
print('NPY data read.')

# read param
global param
param = np.load(runpath+'/param.npy').all()
param['dat_type'] = np.float32

# import functions
exec(open(path+'swm_param.py').read())
exec(open(path+'swm_operators.py').read())
exec(open(path+'swm_output.py').read())
param['output'] = 0

set_grad_mat()
set_interp_mat()
set_lapl_mat()
set_coriolis()

## create ouputfolder
try:
    os.mkdir(runpath+'/plots')
except:
   pass
   
## PLOTTING   
hm_max = np.percentile(abs(hm),98)
levs = np.linspace(-hm_max,hm_max,64)
    
fig,ax = plt.subplots(1,1,figsize=(12,9))

c = ax.contourf(param['x_T']/1e3,param['y_T']/1e3,h2mat(hm),levs,cmap='RdBu_r',extend='both')    
plt.colorbar(c,ax=ax)

xx_u,yy_u = np.meshgrid(param['x_u']/1e3,param['y_u']/1e3)

# reduce the number of plotted quiver arrows in each dim to 60
nquivx = int(param['nx']/min(60,param['nx']))    
nquivy = int(param['ny']/min(60,param['ny']))
qs = np.ogrid[0:param['ny']:nquivy,0:param['nx']-1:nquivx]

ax.quiver(xx_u[qs],yy_u[qs],u2mat(um)[qs],u2mat(Ivu.dot(vm))[qs])

ax.set_ylabel('y [km]')
ax.set_xlabel('x [km]')
ax.set_title(r'Mean state $u,v,\eta$')

plt.tight_layout()
plt.savefig(runpath+'/plots/uvh_mean.png')
plt.close(fig)

# MEAN KINETIC ENERGY
levs = np.linspace(0,np.percentile(mke,98),64)

fig,ax = plt.subplots(1,1,figsize=(12,9))

c1 = ax.contourf(param['x_T']/1e3,param['y_T']/1e3,h2mat(mke),levs,extend='max',cmap='viridis')    
cb = plt.colorbar(c1,ax=ax)
cb.set_label(r'Jm$^{-2}$')

ax.set_ylabel('y [km]')
ax.set_xlabel('x [km]')
ax.set_title(r'Mean kinetic energy')

plt.tight_layout()
plt.savefig(runpath+'/plots/mke_mean.png')
plt.close(fig)

# Mean Potential Energy
levs = np.linspace(0,np.percentile(mpe,98),64)

fig,ax = plt.subplots(1,1,figsize=(12,9))

c1 = ax.contourf(param['x_T']/1e3,param['y_T']/1e3,h2mat(mpe),levs,extend='max',cmap='viridis')    
cb = plt.colorbar(c1,ax=ax)
cb.set_label(r'Jm$^{-2}$')

ax.set_ylabel('y [km]')
ax.set_xlabel('x [km]')
ax.set_title(r'Mean potential energy')

plt.tight_layout()
plt.savefig(runpath+'/plots/mpe_mean.png')
plt.close(fig)

## VARIANCE
# U VARIANCE
levs = np.linspace(0,np.percentile(uvar,98),64)

fig,ax = plt.subplots(1,1,figsize=(12,9))

c1 = ax.contourf(param['x_u']/1e3,param['y_u']/1e3,u2mat(uvar),levs,extend='max',cmap='viridis')    
cb = plt.colorbar(c1,ax=ax)
cb.set_label(r'm$^2$s$^{-2}$')

ax.set_ylabel('y [km]')
ax.set_xlabel('x [km]')
ax.set_title(r'Variance $u$')

plt.tight_layout()
plt.savefig(runpath+'/plots/u_var.png')
plt.close(fig)

# V VARIANCE
levs = np.linspace(0,np.percentile(vvar,98),64)

fig,ax = plt.subplots(1,1,figsize=(12,9))

c1 = ax.contourf(param['x_v']/1e3,param['y_v']/1e3,v2mat(vvar),levs,extend='max',cmap='viridis')    
cb = plt.colorbar(c1,ax=ax)
cb.set_label(r'm$^2$s$^{-2}$')

ax.set_ylabel('y [km]')
ax.set_xlabel('x [km]')
ax.set_title(r'Variance $v$')

plt.tight_layout()
plt.savefig(runpath+'/plots/v_var.png')
plt.close(fig)

# Eddy Potential Energy (proportional to h variance)
levs = np.linspace(0,np.percentile(epe,98),64)

fig,ax = plt.subplots(1,1,figsize=(12,9))

c1 = ax.contourf(param['x_T']/1e3,param['y_T']/1e3,h2mat(epe),levs,extend='max',cmap='viridis')    
cb = plt.colorbar(c1,ax=ax)
cb.set_label(r'Jm$^{-2}$')

ax.set_ylabel('y [km]')
ax.set_xlabel('x [km]')
ax.set_title(r'Eddy potential energy')

plt.tight_layout()
plt.savefig(runpath+'/plots/epe_mean.png')
plt.close(fig)

# EKE
levs = np.linspace(0,np.percentile(eke,98),64)

fig,ax = plt.subplots(1,1,figsize=(12,9))

c1 = ax.contourf(param['x_T']/1e3,param['y_T']/1e3,h2mat(eke),levs,extend='max',cmap='viridis')    
cb = plt.colorbar(c1,ax=ax)
cb.set_label(r'Jm$^{-2}$')

ax.set_ylabel('y [km]')
ax.set_xlabel('x [km]')
ax.set_title(r'Eddy kinetic energy')

plt.tight_layout()
plt.savefig(runpath+'/plots/eke_mean.png')
plt.close(fig)

## RELATIVE VORTICITY MEAN
zm_max = np.percentile(abs(zm),98)
levs = np.linspace(-zm_max,zm_max,64)

fig,ax = plt.subplots(1,1,figsize=(10,9))
 
ax.contourf(param['x_q']/1e3,param['y_q']/1e3,q2mat(zm),levs,cmap='RdBu_r',extend='both')    

ax.set_xlim(0,param['Lx']/1e3)
ax.set_ylim(0,param['Ly']/1e3)

ax.set_ylabel('y [km]')
ax.set_xlabel('x [km]')
ax.set_title(r'Mean relative vorticity')

plt.tight_layout()
plt.savefig(runpath+'/plots/relvort_mean.png')
plt.close(fig)

## SPEED MEAN
levs = np.linspace(0,np.percentile(speedm,99),64)

fig,ax = plt.subplots(1,1,figsize=(12,9))

c1 = ax.contourf(param['x_T']/1e3,param['y_T']/1e3,h2mat(speedm),levs,extend='max',cmap='viridis')    
cb = plt.colorbar(c1,ax=ax)
cb.set_label(r'$|\mathbf{u}|$')

ax.set_ylabel('y [km]')
ax.set_xlabel('x [km]')
ax.set_title('Mean speed')

plt.tight_layout()
plt.savefig(runpath+'/plots/speed_mean.png')
plt.close(fig)

## REYNOLDS
Re = np.log10(Re)
levs = np.linspace(np.percentile(Re,2),np.percentile(Re,98),64)

fig,ax = plt.subplots(1,1,figsize=(12,9))

c1 = ax.contourf(param['x_T']/1e3,param['y_T']/1e3,h2mat(Re),levs,extend='both',cmap='viridis')    
cb = plt.colorbar(c1,ax=ax)
cb.set_label(r'log$_{10}$(Re)')

ax.set_ylabel('y [km]')
ax.set_xlabel('x [km]')
ax.set_title('Mean Reynolds number')

plt.tight_layout()
plt.savefig(runpath+'/plots/Re_mean.png')
plt.close(fig)

## ROSSBY
Ro = np.log10(Ro)
levs = np.linspace(np.percentile(Ro,2),np.percentile(Ro,98),64)

fig,ax = plt.subplots(1,1,figsize=(12,9))

c1 = ax.contourf(param['x_T']/1e3,param['y_T']/1e3,h2mat(Ro),levs,extend='both',cmap='viridis')    
cb = plt.colorbar(c1,ax=ax)
cb.set_label(r'log$_{10}$(Ro)')

ax.set_ylabel('y [km]')
ax.set_xlabel('x [km]')
ax.set_title('Mean Rossby number')

plt.tight_layout()
plt.savefig(runpath+'/plots/Ro_mean.png')
plt.close(fig)

## EKMAN
Ek = np.log10(Ek)
levs = np.linspace(np.percentile(Ek,2),np.percentile(Ek,98),64)

fig,ax = plt.subplots(1,1,figsize=(12,9))

c1 = ax.contourf(param['x_T']/1e3,param['y_T']/1e3,h2mat(Ek),levs,extend='both',cmap='viridis')    
cb = plt.colorbar(c1,ax=ax)
cb.set_label(r'log$_{10}$(Ek)')

ax.set_ylabel('y [km]')
ax.set_xlabel('x [km]')
ax.set_title('Mean Ekman number')

plt.tight_layout()
plt.savefig(runpath+'/plots/Ek_mean.png')
plt.close(fig)

## TIME SERIES OF NONDIM NUMBERS
fig,(ax1,ax2,ax3) = plt.subplots(3,1,sharex=True,figsize=(8,5))

ax1.plot(t,Rem)
ax2.plot(t,Rom)
ax3.plot(t,Ekm)

ax1.set_xlim(t[0],t[-1])
ax1.set_title('Reynolds number')
ax2.set_title('Rossby number')
ax3.set_title('Ekman number')
ax3.set_xlabel('time [days]')

plt.tight_layout()
plt.savefig(runpath+'/plots/ReRoEk_tseries.png')
plt.close(fig)

## ENERGY AND ENSTROPHY TIME SERIES
fig,(ax1,ax2,ax3) = plt.subplots(3,1,sharex=True,figsize=(8,5))

# do not show initial conditions
ax1.plot(t,KEm)
ax2.plot(t,PEm)
ax3.plot(t,Zm)

ax1.set_title('kinetic energy')
ax2.set_title('potential energy')
ax3.set_title('potential enstrophy')

ax3.set_xlabel('time [days]')

plt.tight_layout()
plt.savefig(runpath+'/plots/Energy_tseries.png')
plt.close(fig)