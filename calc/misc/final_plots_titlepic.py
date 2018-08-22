## PRODUCE MEAN PLOTS
from __future__ import print_function
#path = '/home/mkloewer/python/swm/'
path = '/network/aopp/cirrus/pred/kloewer/swm_superhr/'
import os; os.chdir(path) # change working directory
import numpy as np
from scipy import sparse
from scipy.integrate import cumtrapz
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import time as tictoc
from netCDF4 import Dataset
import glob
from cmocean import cm
from matplotlib.colors import BoundaryNorm

## OPTIONS
runfolder = 16
print(('Creating final plots from run %i') % runfolder)

## read data
runpath = path+'run%04i' % runfolder
ncu = Dataset(runpath+'/u.nc')
ncv = Dataset(runpath+'/v.nc')
#nch = Dataset(runpath+'/h.nc')

u = ncu['u'][-1,:,:].flatten()
v = ncv['v'][-1,:,:].flatten()
#h = nch['h'][-1,:,:].flatten()
#time = nch['t'][:][:]
print('netCDF data read.')

# close netcdfs
ncu.close()
ncv.close()
#nch.close()

# read param
global param

#param = dict()
#param['nx'] = 1024
#param['ny'] = 1024

#param['NT'] = param['nx']*param['ny']
#param['Nu'] = (param['nx']-1)*param['ny']
#param['Nv'] = param['nx']*(param['ny']-1)
#param['Nq'] = (param['nx']+1)*(param['ny']+1)

param = np.load(runpath+'/param.npy').all()
#param['t_end'] = time[-1]/(3600.*24.)
param['dat_type'] = np.float32

path2 = '/a/cplxfs3/srv/homes/aopp/kloewer/swm_opt_cntrl/'
path2 = '/home/kloewer/swm_opt_cntrl/'

# import functions
exec(open(path2+'swm_param.py').read())
exec(open(path2+'swm_operators.py').read())
exec(open(path2+'swm_output.py').read())
exec(open(path2+'swm_rhs.py').read())
param['output'] = 0

set_grad_mat()
set_interp_mat()
set_lapl_mat()
#set_coriolis()

## create ouputfolder
try:
    os.mkdir(runpath+'/plots')
except:
   pass
    
## U,V,H final
# h_max = np.percentile(abs(h),99)
# levs = np.linspace(-h_max,h_max,64)
#     
# fig,ax = plt.subplots(1,1,figsize=(12,9))
# 
# c = ax.contourf(param['x_T']/1e3,param['y_T']/1e3,h2mat(h),levs,cmap='RdBu_r',extend='both')    
# plt.colorbar(c,ax=ax)
# 
# xx_u,yy_u = np.meshgrid(param['x_u']/1e3,param['y_u']/1e3)
# 
# #reduce the number of plotted quiver arrows in each dim to 60
# nquivx = int(param['nx']/min(60,param['nx']))    
# nquivy = int(param['ny']/min(60,param['ny']))
# qs = np.ogrid[0:param['ny']:nquivy,0:param['nx']-1:nquivx]
# 
# ax.quiver(xx_u[qs],yy_u[qs],u2mat(u)[qs],u2mat(Ivu.dot(v))[qs])
# 
# ax.set_ylabel('y [km]')
# ax.set_xlabel('x [km]')
# 
# ax.set_title(r'$u,v,\eta$ at $t$ = %.2f days' % param['t_end'])
# 
# plt.tight_layout()
# plt.savefig(runpath+'/plots/uvh_final.png')
# plt.close(fig)

## REL VORTICITY FINAL
z = Gvx.dot(v) - Guy.dot(u)
#z = np.sqrt(ITq.dot((Gux.dot(u) - Gvy.dot(v))**2) + (Guy.dot(u) + Gvx.dot(v))**2)
z = np.sign(z)*abs(z)**(1/2.)

z_max = np.percentile(abs(z),98.5)
levs = np.linspace(-z_max,z_max,128)

fig,ax = plt.subplots(1,1,figsize=(9,9))

q = ax.contourf(param['x_q']/1e3,param['y_q']/1e3,q2mat(z),levs,cmap=cm.balance,extend='both')    

#plt.colorbar(q,ax=ax)

ax.set_xlim(0,param['Lx']/1e3)
ax.set_ylim(0,param['Ly']/1e3)

ax.set_xticks([])
ax.set_yticks([])

#ax.set_ylabel('y [km]')
#ax.set_xlabel('x [km]')
#ax.set_title(r'Relative vorticity at $t$ = %.2f days' % param['t_end'])

plt.tight_layout()
plt.savefig(runpath+'/plots/relvort_final.png',dpi=300)
plt.close(fig)

# 
# ## LENGTH SCALE FINAL
# z = Gvx.dot(v) - Guy.dot(u)
# KE = .5*(IuT.dot(u**2) + IvT.dot(v**2))
# L = np.sqrt(np.sqrt(KE/IqT.dot(z**2)))
# 
# L_max = np.percentile(L,98)
# levs = np.linspace(0,L_max,64)
# 
# fig,ax = plt.subplots(1,1,figsize=(10,9))
# 
# q = ax.contourf(param['x_T']/1e3,param['y_T']/1e3,h2mat(L),levs,cmap='viridis',extend='max')    
# 
# plt.colorbar(q,ax=ax)
# ax.set_xlim(0,param['Lx']/1e3)
# ax.set_ylim(0,param['Ly']/1e3)
# 
# ax.set_ylabel('y [km]')
# ax.set_xlabel('x [km]')
# ax.set_title(r'Length scale at $t$ = %.2f days' % param['t_end'])
# 
# plt.tight_layout()
# plt.savefig(runpath+'/plots/length_final.png')
# plt.savefig(runpath+'/plots/length_final.png')
# plt.close(fig)
#     
# 
# ## SPEED FINAL
# speed = np.sqrt(IuT.dot(u**2) + IvT.dot(v**2))
# 
# levs = np.linspace(0,np.percentile(speed,99),64)
# 
# fig,ax = plt.subplots(1,1,figsize=(12,9))
# 
# c1 = ax.contourf(param['x_T']/1e3,param['y_T']/1e3,h2mat(speed),levs,extend='max',cmap='viridis')    
# cb = plt.colorbar(c1,ax=ax)
# cb.set_label(r'$|\mathbf{u}|$')
# 
# ax.set_ylabel('y [km]')
# ax.set_xlabel('x [km]')
# ax.set_title(r'Speed [m/s], $t$=%.2f days' % param['t_end'])
# 
# plt.tight_layout()
# plt.savefig(runpath+'/plots/speed_final.png')
# plt.close(fig)
#     
# ## REYNOLDS, ROSSBY AND EKMAN FINAL
# u_T = IuT.dot(u)    # u on T-grid
# v_T = IvT.dot(v)    # v on T-grid
# 
# #advective term
# adv_u = u_T*Gux.dot(u) + v_T*IqT.dot(Guy.dot(u))
# adv_v = u_T*IqT.dot(Gvx.dot(v)) + v_T*Gvy.dot(v)
# adv_term = np.sqrt(adv_u**2 + adv_v**2)
# 
# #coriolis term
# cor_term = f_T*np.sqrt(IuT.dot(u**2) + IvT.dot(v**2))
# 
# #diffusive term
# diff_u = param['nu_B']*LLu.dot(u)
# diff_v = param['nu_B']*LLv.dot(v)
# diff_term = np.sqrt(IuT.dot(diff_u**2) + IvT.dot(diff_v**2))
# 
# Re = adv_term / diff_term
# Ek = diff_term / cor_term
# Ro = adv_term / cor_term
# """
# diff_u2,diff_v2 = mixing(u,v,h+500.,ITq.dot(h+500.),ITu.dot(h+500.),ITv.dot(h+500.))
# diff_u2 = IuT.dot(diff_u2)
# diff_v2 = IvT.dot(diff_v2)
# Re2 = np.sqrt(adv_u**2/diff_u2**2 + adv_v**2/diff_v2**2)
# Re2m = np.logical_or((diff_u2 == 0),(diff_v2 == 0))
# Re2 = np.ma.masked_array(np.log10(Re2),mask=Re2m)
# """
# # REYNOLDS
# fig,ax = plt.subplots(1,1,figsize=(12,9))
# 
# Re = np.log10(Re) # use the logarithm instead
# Re_min = np.percentile(Re,3)
# Re_max = np.percentile(Re,99)
# levs = np.linspace(Re_min,Re_max,64)
# 
# c = ax.contourf(param['x_T']/1e3,param['y_T']/1e3,h2mat(Re),levs,extend='both',cmap='viridis') 
# cb = plt.colorbar(c,ax=ax)
# cb.set_label(r'log$_{10}$(Re)')
# 
# ax.set_ylabel('y [km]')
# ax.set_xlabel('x [km]')
# ax.set_title('Reynolds number, $t$=%.2f days' % param['t_end'])
# 
# plt.tight_layout()
# plt.savefig(runpath+'/plots/Re_final.png')
# plt.close(fig)
# #plt.show()
# """
# # REYNOLDS 2
# fig,ax = plt.subplots(1,1,figsize=(12,9))
# 
# #Re = np.log10(Re) # use the logarithm instead
# #Re_min = np.percentile(Re,2)
# #Re_max = np.percentile(Re,98)
# #levs = np.linspace(Re_min,Re_max,64)
# 
# c = ax.contourf(param['x_T']/1e3,param['y_T']/1e3,h2mat(Re2),levs,extend='both',cmap='viridis') 
# cb = plt.colorbar(c,ax=ax)
# cb.set_label(r'log$_{10}$(Re)')
# 
# ax.set_ylabel('y [km]')
# ax.set_xlabel('x [km]')
# ax.set_title('Reynolds number, $t$=%.2f days' % param['t_end'])
# 
# plt.tight_layout()
# plt.savefig(runpath+'/plots/Re_final2.png')
# plt.close(fig)
# """
# ##
# 
# # ROSSBY PLOT
# fig,ax = plt.subplots(1,1,figsize=(12,9))
# 
# Ro = np.log10(Ro)
# Ro_min = np.percentile(Ro,2)
# Ro_max = np.percentile(Ro,98)
# levs = np.linspace(Ro_min,Ro_max,64)
# 
# c = ax.contourf(param['x_T']/1e3,param['y_T']/1e3,h2mat(Ro),levs,extend='both',cmap='viridis') 
# cb = plt.colorbar(c,ax=ax)
# cb.set_label(r'log$_{10}$(Ro)')
# 
# ax.set_ylabel('y [km]')
# ax.set_xlabel('x [km]')
# ax.set_title('Rossby number, $t$=%.2f days' % param['t_end'])
# 
# plt.tight_layout()
# plt.savefig(runpath+'/plots/Ro_final.png')
# plt.close(fig)
#     
# # EKMAN PLOT
# fig,ax = plt.subplots(1,1,figsize=(12,9))
# 
# Ek = np.log10(Ek)
# Ek_min = np.percentile(Ek,2)
# Ek_max = np.percentile(Ek,98)
# levs = np.linspace(Ek_min,Ek_max,64)
# 
# c = ax.contourf(param['x_T']/1e3,param['y_T']/1e3,h2mat(Ek),levs,extend='both',cmap='viridis') 
# cb = plt.colorbar(c,ax=ax)
# cb.set_label(r'log$_{10}$(Ek)')
# 
# ax.set_ylabel('y [km]')
# ax.set_xlabel('x [km]')
# ax.set_title('Ekman number, $t$=%.2f days' % param['t_end'])
# 
# plt.tight_layout()
# plt.savefig(runpath+'/plots/Ek_final.png')
