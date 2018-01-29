## PRODUCE MEAN CALCULATIONS AND EXPORT AS .NPY
from __future__ import print_function

# path
import os
funpath = '/network/home/aopp/kloewer/git/swm/'
path = '/network/aopp/cirrus/pred/kloewer/swm_back_ronew/'
os.chdir(path)                              # change working directory

import numpy as np
from scipy import sparse

# OPTIONS
runfolder = [5]
print('Calculating timeseries from run ' + str(runfolder))

## read data
for r,i in zip(runfolder,range(len(runfolder))):
    runpath = path+'run%04i' % r

    if i == 0:
        #skip = 5*365
        skip = 0
        u = np.load(runpath+'/u_sub.npy')[skip:,...]
        v = np.load(runpath+'/v_sub.npy')[skip:,...]
        eta = np.load(runpath+'/eta_sub.npy')[skip:,...]
        t = np.load(runpath+'/t_sub.npy')[skip:,...]
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
exec(open(funpath+'swm_param.py').read())
exec(open(funpath+'swm_operators.py').read())

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
# thickness-weighted averaging for u,v (Aiki, 2016)
#etam = eta.mean(axis=1)             # temporal average
#um = (ITu.dot(eta + param['H'])*u).mean(axis=1) / ITu.dot(etam + param['H'])
#vm = (ITv.dot(eta + param['H'])*v).mean(axis=1) / ITv.dot(etam + param['H'])
#print('u,v,eta mean done.')

# TOTAL KINETIC ENERGY
#ke = .5*param['rho']*((eta+param['H'])*(IuT.dot(u**2) + IvT.dot(v**2))).mean(axis=1)
# TOTAL (PERTURBATION) POTENTIAL ENERGY
#pe = .5*param['rho']*param['g']*(eta**2).mean(axis=1)
#print('KE, PE done.')

# MEAN KINETIC ENERGY
#mke = .5*param['rho']*(etam+param['H'])*(IuT.dot(um**2) + IvT.dot(vm**2))
# Mean Potential Energy
#mpe = .5*param['g']*param['rho']*etam**2
#print('MKE, MPE done.')

# Eddy potential energy
#epe = .5*param['rho']*param['g']*eta.var(axis=1)

# Eddy kinetic energy
#uprime = ((u.T - um).T)
#vprime = ((v.T - vm).T)
#print('Prime anomalies done.')

#eke = .5*param['rho']*((eta+param['H'])*(IuT.dot(uprime**2) + IvT.dot(vprime**2))).mean(axis=1)
#print('EPE, EKE done.')
#del uprime,vprime

# SPEED MEAN
#speedm = np.sqrt(IuT.dot(u**2) + IvT.dot(v**2)).mean(axis=1)
#print('Rel Vort, Speed done.')

## ENERGY TIME SERIES
PEm = (.5*param['g']*param['rho']*eta**2).mean(axis=0)
print('Potential Energy done.')

KEm = .5*param['rho']*((eta+param['H'])*(IuT.dot(u**2) + IvT.dot(v**2))).mean(axis=0)
print('Kinetic Energy done.')

## STORING
dic = dict()
#all_var2export = ['um','vm','etam','ke','pe','mke','eke','mpe','epe','speedm']
all_var2export = ['t','PEm','KEm']

for v in all_var2export:
    exec('dic[v] ='+v)

np.save(runpath+'/analysis/mean_timeseries.npy',dic)
print('Everything stored.')
