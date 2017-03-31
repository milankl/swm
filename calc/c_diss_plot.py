from __future__ import print_function
path = '/home/mkloewer/python/swm/'
import os; os.chdir(path) # change working directory
import numpy as np
from scipy import sparse
import time as tictoc
from netCDF4 import Dataset
import glob
import matplotlib.pyplot as plt

# OPTIONS
runfolder = [2,3]

## read data
for r,i in zip(runfolder,range(len(runfolder))):
    runpath = path+'data/run%04i' % r
    
    if i == 0:
        u = np.load(runpath+'/u_sub.npy')
        v = np.load(runpath+'/v_sub.npy')
        h = np.load(runpath+'/h_sub.npy')
        time = np.load(runpath+'/t_sub.npy')
        print('run %i read.' % r)

    else:
        u = np.concatenate((u,np.load(runpath+'/u_sub.npy')))
        v = np.concatenate((v,np.load(runpath+'/v_sub.npy')))
        h = np.concatenate((h,np.load(runpath+'/h_sub.npy')))
        time = np.hstack((time,np.load(runpath+'/t_sub.npy')))
        print('run %i read.' % r)

t = time / 3600. / 24.  # in days
## read param
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
    os.mkdir(runpath+'/analysis')
except:
   pass
   
## reshape u,v
u = u.reshape((tlen,param['Nu'])).T
v = v.reshape((tlen,param['Nv'])).T
h = h.reshape((tlen,param['NT'])).T
print('Reshape done.')

##
dudx = Gux.dot(u)
dudy = Guy.dot(u)
dvdx = Gvx.dot(v)
dvdy = Gvy.dot(v)

n = 2

D = np.sqrt((dudx - dvdy)**2 + IqT.dot((dudy + dvdx)**2))
Ro = (D.T/f_T)
Rom = Ro.mean(axis=0)
c = (1/(1+Ro)**n).mean(axis=0)

# REYNOLDS, ROSSBY, EKMAN NUMBER MEAN
u_T = IuT.dot(u)
v_T = IvT.dot(v)
print('u,v interpolation done.')

#advective term
adv_u = u_T*Gux.dot(u) + v_T*IqT.dot(Guy.dot(u))
adv_v = u_T*IqT.dot(Gvx.dot(v)) + v_T*Gvy.dot(v)
del u_T,v_T
adv_term = np.sqrt(adv_u**2 + adv_v**2)
del adv_u, adv_v
print('Advection term done.')

#coriolis term
cor_term = (f_T*np.sqrt(IuT.dot(u**2) + IvT.dot(v**2)).T).T
print('Coriolis term done.')

Ro2 = adv_term / cor_term
c2 = (1/(1+Ro2)**n).mean(axis=1)
Ro2m = Ro2.mean(axis=1)

##
levs1 = np.linspace(0,.2,21)
levs2 = np.linspace(0.5,1,21)

fig,axs = plt.subplots(2,3,sharex=True,sharey=True,figsize=(9,5.5))
plt.tight_layout(rect=[-.02,-.03,1.12,.97],w_pad=0.1)

axs[0,0].contourf(param['x_T'],param['y_T'],h2mat(Ro2m),levs1)
axs[0,1].contourf(param['x_T'],param['y_T'],h2mat(Rom),levs1,extend='max')
m1 = axs[0,2].contourf(param['x_T'],param['y_T'],h2mat(Ro[-1,:]),levs1,extend='max')
plt.colorbar(m1,ax=(axs[0,0],axs[0,1],axs[0,2]),ticks=np.arange(0,.22,.04))

axs[1,0].contourf(param['x_T'],param['y_T'],h2mat(c2),levs2)
m21 = axs[1,0].contour(param['x_T'],param['y_T'],h2mat(c2),[0.8],linewidths=0.7)
axs[1,1].contourf(param['x_T'],param['y_T'],h2mat(c),levs2)
m2 = axs[1,2].contourf(param['x_T'],param['y_T'],h2mat(1/(1+Ro[-1,:])**n),levs2,extend='min')
axs[1,2].contour(param['x_T'],param['y_T'],h2mat(1/(1+Ro[-1,:])**n),[0.8],linewidths=0.7)
m22 = axs[1,1].contour(param['x_T'],param['y_T'],h2mat(c),[0.8],linewidths=0.7)
plt.colorbar(m2,ax=(axs[1,0],axs[1,1],axs[1,2]),ticks=np.arange(0.5,1.05,.05))
plt.clabel(m22, inline=1, fontsize=5,fmt='%.1f')
plt.clabel(m21, inline=1, fontsize=5,fmt='%.1f')

axs[0,0].set_xticks([])
axs[0,0].set_yticks([])

axs[0,0].set_title(r'$\overline{R_o} = \overline{\frac{|(\mathbf{u} \cdot \nabla)\mathbf{u}|}{|f\mathbf{u}|}}$')
axs[0,1].set_title(r'$\overline{R_o^*} = \overline{\frac{|D|}{f}}$')
axs[0,2].set_title(r'snapshot: $R_o^*$')

axs[1,0].set_title(r'$(1+\overline{R_o})^{-2}$')
axs[1,1].set_title(r'$(1+\overline{R_o}^*)^{-2}$')
axs[1,2].set_title(r'$(1+R_o^*)^{-2}$')

axs[0,0].set_ylabel('y')
axs[1,0].set_ylabel('y')
axs[1,0].set_xlabel('x')
axs[1,1].set_xlabel('x')

plt.savefig(path+'compare/Ro_scaling.png',dpi=150)
plt.close(fig)
#plt.show()


