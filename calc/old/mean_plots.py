## PRODUCE MEAN PLOTS
from __future__ import print_function
path = '/home/mkloewer/python/swm/'
import os; os.chdir(path) # change working directory
import numpy as np
from scipy import sparse
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time as tictoc
from netCDF4 import Dataset
import glob
from matplotlib.colors import BoundaryNorm,LogNorm

# OPTIONS
runfolder = 5
print(('Creating mean plots from run %i') % runfolder)

## read data
runpath = path+'stoch/data/run%04i' % runfolder

ncu = Dataset(runpath+'/u.nc')
u = ncu['u'][:][:,:,:]
ncu.close()
print('u read.')

ncv = Dataset(runpath+'/v.nc')
v = ncv['v'][:][:,:,:]
ncv.close()
print('v read.')

nch = Dataset(runpath+'/h.nc')
h = nch['h'][:][:,:,:]
time = nch['t'][:]   # in seconds
t = time / 3600. / 24.  # in days
nch.close()
print('netCDF data read.')

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

tlen = len(time)
## create ouputfolder
try:
    os.mkdir(runpath+'/plots')
except:
   pass
   
## reshape u,v
u = u.reshape((tlen,param['Nu'])).T
v = v.reshape((tlen,param['Nv'])).T
h = h.reshape((tlen,param['NT'])).T
print('Reshape done.')
    
## U,V,H mean
um = u.mean(axis=1)  # temporal averages
vm = v.mean(axis=1)
hm = h.mean(axis=1)

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

print('uvh mean done.')
del xx_u,yy_u

# MEAN KINETIC ENERGY
mke = .5*param['rho']*param['H']*(IuT.dot(um**2) + IvT.dot(vm**2))
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
mpe = .5*param['g']*param['rho']*hm**2
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
uvar = u.var(axis=1)  # temporal averages
vvar = v.var(axis=1)
epe = h.var(axis=1)*param['g']*param['rho']/2.   # eddy potential energy proportional to hvar

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
eke = .5*param['rho']*param['H']*(IuT.dot(uvar) + IvT.dot(vvar))
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

print('Variance done.')
del uvar,vvar

## RELATIVE VORTICITY MEAN
tlen = u.shape[0]   # length of time dimension
z = Gvx.dot(v) - Guy.dot(u)
print('Vorticity calculated')
zm = z.mean(axis=1)  # temporal average
zm = np.sign(zm)*abs(zm)**(1/2.) # plot actually the sqrt of rel vort

zm_max = np.percentile(abs(zm),98)
levs = np.linspace(-zm_max,zm_max,64)

fig,ax = plt.subplots(1,1,figsize=(10,9))

#ax.pcolormesh(param['x_q']/1e3,param['y_q']/1e3,q2mat(zm),cmap='RdBu_r',norm=BoundaryNorm(levs,ncolors=256))    
ax.contourf(param['x_q']/1e3,param['y_q']/1e3,q2mat(zm),levs,cmap='RdBu_r',extend='both')    

ax.set_xlim(0,param['Lx']/1e3)
ax.set_ylim(0,param['Ly']/1e3)

ax.set_ylabel('y [km]')
ax.set_xlabel('x [km]')
ax.set_title(r'Mean relative vorticity')

plt.tight_layout()
plt.savefig(runpath+'/plots/relvort_mean.png')
plt.close(fig)
print('Vorticiy done.')
del z,zm

## SPEED MEAN
speedm = np.sqrt(IuT.dot(u**2) + IvT.dot(v**2)).mean(axis=1)
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

print('Speed done.')
del speedm

## REYNOLDS, ROSSBY, EKMAN NUMBER MEAN
u_T = IuT.dot(u)
v_T = IvT.dot(v)
print('u,v interpolation done.')

#advective term
adv_u = u_T*Gux.dot(u) + v_T*IqT.dot(Guy.dot(u))
adv_v = u_T*IqT.dot(Gvx.dot(v)) + v_T*Gvy.dot(v)
del u_T,v_T
adv_term = np.sqrt(adv_u**2 + adv_v**2)
del adv_u, adv_v
adv_termm = adv_term.mean(axis=0)   # spatial average
adv_term = adv_term.mean(axis=1)    # temporal average
print('Advection term done.')

#coriolis term
cor_term = (f_T*np.sqrt(IuT.dot(u**2) + IvT.dot(v**2)).T).T
cor_termm = cor_term.mean(axis=0)   # spatial average
cor_term = cor_term.mean(axis=1)    # temporal average
print('Coriolis term done.')

#diffusive term
diff_u = param['nu_B']*LLu.dot(u)
diff_v = param['nu_B']*LLv.dot(v)
diff_term = np.sqrt(IuT.dot(diff_u**2) + IvT.dot(diff_v**2))
del diff_u,diff_v
diff_termm = diff_term.mean(axis=0)    # spatial average
diff_term = diff_term.mean(axis=1)     # temporal average
print('Diffusion term done.')

# temporal averages
Re = np.log10(adv_term / diff_term)
Ek = np.log10(diff_term / cor_term)
Ro = np.log10(adv_term / cor_term)

# spatial averages
Rem = adv_termm / diff_termm
Ekm = diff_termm / cor_termm
Rom = adv_termm / cor_termm
print('Nondim numbers done.')

## REYNOLDS
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

print('Nondim numbers plots done')

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
PEm = (.5*param['g']*param['rho']*h**2).mean(axis=0)
print('Potential Energy done.')

KEm = (.5*param['rho']*param['H']*(IuT.dot(u**2) + IvT.dot(v**2))).mean(axis=0)
print('Kinetic Energy done.')

PV = (f_q + (Gvx.dot(v) - Guy.dot(u)).T).T / ITq.dot(h+param['H'])
Zm = (PV**2*ITq.dot(h+param['H'])).mean(axis=0)
del PV
print('Enstrophy done.')

fig,(ax1,ax2,ax3) = plt.subplots(3,1,sharex=True,figsize=(8,5))

# do not show initial conditions
ax1.plot(t,KEm)
ax2.plot(t,PEm)
ax3.plot(t,Zm)

ax1.set_title(r'kinetic energy / KE$_1$')
ax2.set_title(r'potential energy / PE$_1$')
ax3.set_title(r'potential enstrophy / Z$_1$')

ax3.set_xlabel('time [days]')

plt.tight_layout()
plt.savefig(runpath+'/plots/Energy_tseries.png')
plt.close(fig)