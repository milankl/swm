## POWER INPUT EXIT COMPUTION AND PLOTTING
from __future__ import print_function
path = '/home/mkloewer/python/swm/'
import os; os.chdir(path) # change working directory
import numpy as np
from scipy import sparse
import time as tictoc
from netCDF4 import Dataset
import glob
import matplotlib.pyplot as plt

# import functions
exec(open(path+'swm_param.py').read())
exec(open(path+'swm_operators.py').read())
exec(open(path+'swm_output.py').read())
exec(open(path+'swm_rhs.py').read())

# OPTIONS
runfolder = [1,5]
print('Calculating input energy from run ' + str(runfolder))

## read data
runpath1 = path+'data/run%04i' % runfolder[0]
u1 = np.load(runpath1+'/u_sub.npy')[:1800,...]   # read only first 4-5 years
v1 = np.load(runpath1+'/v_sub.npy')[:1800,...]
t1 = np.load(runpath1+'/t_sub.npy')[:1800,...] 
param1 = np.load(runpath1+'/param.npy').all()
D1 = np.load(runpath1+'/analysis/mean.npy').all()
param1['output'] = 0
print('Data from run %i read.' % runfolder[0])

runpath2 = path+'data/run%04i' % runfolder[1]
u2 = np.load(runpath2+'/u_sub.npy')[:1800,...]
v2 = np.load(runpath2+'/v_sub.npy')[:1800,...]
t2 = np.load(runpath2+'/t_sub.npy')[:1800,...]
param2 = np.load(runpath2+'/param.npy').all()
D2 = np.load(runpath2+'/analysis/mean.npy').all()
param2['output'] = 0
print('Data from run %i read.' % runfolder[1])

dt1 = t1[1]-t1[0]   # in seconds
dt2 = t2[1]-t2[0]

t1 = t1 / 3600. / 24. / 365. # in years
t2 = t2 / 3600. / 24. / 365. # in years
## reshape
u1 = u1.reshape((len(t1),param1['Nu'])).T
u2 = u2.reshape((len(t2),param2['Nu'])).T

v1 = v1.reshape((len(t1),param1['Nv'])).T
v2 = v2.reshape((len(t2),param2['Nv'])).T
print('Reshape done.')

## read param
global param
param = param1 # param1 is now the active param in the set_?() functions
set_forcing()
set_grad_mat()
set_lapl_mat()
In1 = ((u1.T*Fx).T*param['rho']*param['H']).mean(axis=0)
In1i = np.cumsum((u1.T*Fx).T*param['rho']*param['H']*dt1,axis=1).mean(axis=0)
print('Input energy 1 done.')
Ex1u = (param1['B']*(u1*LLu.dot(u1))*param['rho']*param['H']).mean(axis=0)
Ex1v = (param1['B']*(v1*LLv.dot(v1))*param['rho']*param['H']).mean(axis=0)
Ex1ui = np.cumsum(Ex1u*dt1,axis=0)
Ex1vi = np.cumsum(Ex1v*dt1,axis=0)    
print('Exit energy 1 done.')

param = param2  # param2 is now the active param in the set_?() functions
set_forcing()
set_grad_mat()
set_lapl_mat()
In2 = ((u2.T*Fx).T*param['rho']*param['H']).mean(axis=0)
In2i = np.cumsum((u2.T*Fx).T*param['rho']*param['H']*dt2,axis=1).mean(axis=0)
print('Input energy 2 done.')
Ex2u = (param2['B']*(u2*LLu.dot(u2))*param['rho']*param['H']).mean(axis=0)
Ex2v = (param2['B']*(v2*LLv.dot(v2))*param['rho']*param['H']).mean(axis=0)
Ex2ui = np.cumsum(Ex2u*dt2,axis=0)
Ex2vi = np.cumsum(Ex2v*dt2,axis=0)
print('Exit energy 2 done.')

#clear workspace
del u2,v2

## plotting
fig = plt.figure(figsize=(14,8))

ax1 = plt.subplot(221)
ax2 = plt.subplot(223)
ax3 = plt.subplot(122)

ax1.plot(t1,In1,label=r'Low resolution, $\Delta x = 30$km')
ax1.plot(t2,In2,label=r'High resolution, $\Delta x = 7.5$km')
ax1.plot([0,2],[0,0],'k')
ax1.legend(loc=2,fontsize=11,ncol=2)
ax1.set_title('Power input by forcing')
ax1.set_xlim(0,2)
ax1.set_ylim(-0.02,0.04)
ax1.set_ylabel(r'Power [Wm$^{-2}$]')

ax2.plot(t1,Ex1u,'b')
ax2.plot(t2,Ex2u,'g')
ax2.plot(t1,Ex1v,'b',lw=2)
ax2.plot(t2,Ex2v,'g',lw=2)
ax2.set_title('Dissipation power by mixing')
ax2.plot(0,0,'grey',label=r'$u$-component')
ax2.plot(0,0,'grey',lw=2,label=r'$v$-component')
ax2.legend(loc=3,ncol=2,fontsize=11)
ax2.set_xlim(0,2)
ax2.set_ylabel(r'Power [Wm$^{-2}$]')
ax2.set_xlabel(r'time [years]')
ax2.set_ylim(-.01,0)

ax3.plot(t1,In1i*1e-3,'b')
ax3.plot(t2,In2i*1e-3,'g')
ax3.plot(t1,(Ex1ui+Ex1vi)*1e-3,'b--')
ax3.plot(t1,(Ex2ui+Ex2vi)*1e-3,'g--')
ax3.plot(t1,(In1i+Ex1ui+Ex1vi)*1e-3,'b',lw=2)
ax3.plot(t2,(In2i+Ex2ui+Ex2vi)*1e-3,'g',lw=2)

ax3.plot([0,4],[0,0],'k')
#only for legend
ax3.plot(0,0,'b',label=r'Low resolution, $\Delta x = 30$km')
ax3.plot(0,0,'g',label=r'High resolution, $\Delta x = 7.5$km')
ax3.plot(0,0,'grey',label='Energy input by forcing')
ax3.plot(0,0,'grey',ls='--',label='Energy dissipation by mixing')
ax3.plot(0,0,'grey',lw=2,label='Net energy input by forcing and mixing')
ax3.legend(loc=3,fontsize=11)

ax3.set_title('Energy input and dissipation')
ax3.set_xlim(0,4)
ax3.set_ylim(-1e3,1e3)
ax3.set_xlabel(r'time [years]')
ax3.set_ylabel(r'Energy [kJm$^{-2}$]')

plt.tight_layout()
plt.savefig('compare/Input_exit_energy.png')
plt.close(fig)
