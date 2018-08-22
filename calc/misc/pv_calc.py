## PRODUCE MEAN CALCULATIONS AND EXPORT AS .NPY
# FOR POTENTIAL VORTICITY
from __future__ import print_function
path = '/home/mkloewer/python/swm/'
import os; os.chdir(path) # change working directory
import numpy as np
from scipy import sparse
import time as tictoc
from netCDF4 import Dataset
import glob

# OPTIONS
runfolder = [3]
print('Calculating means from run ' + str(runfolder))

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

h_q = ITq.dot(h+param['H'])

## POTENTIAL VORTICITY MEAN
Z = (Gvx.dot(v) - Guy.dot(u))
PVm = ((f_q + Z.T).T/h_q).mean(axis=1)

## LENGTH SCALE MEAN
KE = .5*(IuT.dot((u.T - u.mean(axis=1)).T**2) + IvT.dot((v.T - v.mean(axis=1)).T**2))
Lm = np.sqrt(KE/IqT.dot(Z**2))

## STORING
dic = dict()
all_var2export = ['PVm','Lm']

for var in all_var2export:
    exec('dic[var] ='+var)
    
np.save(runpath+'/analysis/PVm_Lm.npy',dic)
print('Everything stored.')