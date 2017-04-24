## PRODUCE MEAN CALCULATIONS AND EXPORT AS .NPY
from __future__ import print_function

# path
import os
path = os.path.dirname(os.getcwd()) + '/'   # on level above
os.chdir(path)                              # change working directory

import numpy as np
from scipy import sparse

# OPTIONS
runfolder = [2]
print('Calculating means from run ' + str(runfolder))

## read data
for r,i in zip(runfolder,range(len(runfolder))):
    runpath = path+'data/run%04i' % r
    
    if i == 0:
        u = np.load(runpath+'/u_sub.npy')
        v = np.load(runpath+'/v_sub.npy')
        eta = np.load(runpath+'/eta_sub.npy')
        t = np.load(runpath+'/t_sub.npy')
        print('run %i read.' % r)

    else:
        # use [1:,...] to not have one time step double 
        u = np.concatenate((u,np.load(runpath+'/u_sub.npy')))[1:,...]
        v = np.concatenate((v,np.load(runpath+'/v_sub.npy')))[1:,...]
        eta = np.concatenate((eta,np.load(runpath+'/eta_sub.npy')))[1:,...]
        t = np.hstack((t,np.load(runpath+'/t_sub.npy')))[1:,...]
        print('run %i read.' % r)

## read param
global param
param = np.load(runpath+'/param.npy').all()

# import functions
exec(open(path+'swm_param.py').read())
exec(open(path+'swm_operators.py').read())

set_grad_mat()
set_interp_mat()
set_lapl_mat()
set_coriolis()

tlen = len(t)
## create ouputfolder
try:
    os.mkdir(runpath+'/analysis')
except:
   pass
   
## reshape u,v
u = u.reshape((tlen,param['Nu'])).T
v = v.reshape((tlen,param['Nv'])).T
eta = eta.reshape((tlen,param['NT'])).T
print('Reshape done.')


## U,V,H mean
um = u.mean(axis=1)  # temporal averages
vm = v.mean(axis=1)
etam = eta.mean(axis=1)
print('u,v,eta mean done.')

# MEAN KINETIC ENERGY
mke = .5*param['rho']*param['H']*(IuT.dot(um**2) + IvT.dot(vm**2))
# Mean Potential Energy
mpe = .5*param['g']*param['rho']*etam**2
print('MKE, MPE done.')

# VARIANCE
uvar = u.var(axis=1)
vvar = v.var(axis=1)

epe = eta.var(axis=1)*param['g']*param['rho']/2.   # eddy potential energy proportional to hvar
# Eddy kinetic energy
eke = .5*param['rho']*param['H']*(IuT.dot(uvar) + IvT.dot(vvar))
print('Variances, EPE, EKE done.')

# SPEED MEAN
speedm = np.sqrt(IuT.dot(u**2) + IvT.dot(v**2)).mean(axis=1)
print('Rel Vort, Speed done.')

## ENERGY TIME SERIES
PEm = (.5*param['g']*param['rho']*eta**2).mean(axis=0)
print('Potential Energy done.')

KEm = (.5*param['rho']*param['H']*(IuT.dot(u**2) + IvT.dot(v**2))).mean(axis=0)
print('Kinetic Energy done.')

## STORING
dic = dict()
all_var2export = ['um','vm','hm','mke','eke','mpe','epe','speedm']
all_var2export += ['t','PEm','KEm']

for v in all_var2export:
    exec('dic[v] ='+v)
    
np.save(runpath+'/analysis/mean.npy',dic)
print('Everything stored.')

