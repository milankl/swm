## PRODUCE SNAPSHOTS OF LAST TIME STEP
from __future__ import print_function

# path
import os
path = os.path.dirname(os.getcwd()) + '/'   # on level above
os.chdir(path)                              # change working directory

import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt
from cmocean import cm

from netCDF4 import Dataset

## OPTIONS
runfolder = 0
print(('Creating final plots from run %i') % runfolder)

## read data
runpath = path+'data/run%04i' % runfolder
ncu = Dataset(runpath+'/u.nc')
ncv = Dataset(runpath+'/v.nc')
nceta = Dataset(runpath+'/eta.nc')

u = ncu['u'][-1,:,:].flatten()
v = ncv['v'][-1,:,:].flatten()
eta = nceta['eta'][-1,:,:].flatten()
time = nceta['t'][:][:]
print('netCDF data read.')

# close netcdfs
ncu.close()
ncv.close()
nceta.close()

# read param
global param
param = np.load(runpath+'/param.npy').all()

# import functions
exec(open(path+'swm_operators.py').read())

set_grad_mat()
set_interp_mat()
set_lapl_mat()
set_coriolis()

t_end = time[-1]/3600./24. 

## create ouputfolder
try:
    os.mkdir(runpath+'/plots')
except:
   pass
    
## U,V,H final
eta_max = np.percentile(abs(eta),99)
levs = np.linspace(-eta_max,eta_max,64)
    
fig,ax = plt.subplots(1,1,figsize=(12,9))

c = ax.contourf(param['x_T']/1e3,param['y_T']/1e3,h2mat(eta),levs,cmap=cm.balance,extend='both')    
plt.colorbar(c,ax=ax)

xx_u,yy_u = np.meshgrid(param['x_u']/1e3,param['y_u']/1e3)

#reduce the number of plotted quiver arrows in each dim to 60
nquivx = int(param['nx']/min(60,param['nx']))    
nquivy = int(param['ny']/min(60,param['ny']))
qs = np.ogrid[0:param['ny']:nquivy,0:param['nx']-1:nquivx]

ax.quiver(xx_u[qs],yy_u[qs],u2mat(u)[qs],u2mat(Ivu.dot(v))[qs])

ax.set_ylabel('y [km]')
ax.set_xlabel('x [km]')

ax.set_title(r'$u,v,\eta$ at $t$ = %i days' % t_end)

plt.tight_layout()
plt.savefig(runpath+'/plots/uvh_final.png')
plt.close(fig)

## REL VORTICITY FINAL
#z = Gvx.dot(v) - Guy.dot(u)
z = ITq.dot(LT.dot(eta))*param['g']/f_q         # based on laplace(eta)
z = np.sign(z)*abs(z)**0.7                      # non-linear scaling for visualization

z_max = np.percentile(abs(z),98)
levs = np.linspace(-z_max,z_max,64)

fig,ax = plt.subplots(1,1,figsize=(9,9))

ax.contourf(param['x_q']/1e3,param['y_q']/1e3,q2mat(z),levs,cmap=cm.balance,extend='both')    
#plt.colorbar(q,ax=ax)

ax.set_xlim(0,param['Lx']/1e3)
ax.set_ylim(0,param['Ly']/1e3)


ax.set_ylabel('y [km]')
ax.set_xlabel('x [km]')
ax.set_title(r'Relative vorticity at $t$ = %i days' % t_end)

plt.tight_layout()
plt.savefig(runpath+'/plots/relvort_final.png')
plt.close(fig)