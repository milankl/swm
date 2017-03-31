## HISTOGRAM COMPUTATIONS FOR REYNOLDS AND ROSSBY NUMBERS
from __future__ import print_function
path = '/home/mkloewer/python/swm/'
import os; os.chdir(path) # change working directory
import numpy as np
from scipy import sparse
import time as tictoc
from netCDF4 import Dataset
import glob

# OPTIONS
runfolder = [2,3]
print('Calculating histograms from run ' + str(runfolder))

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
    
#TODO CORRECT CALCULATION FOR EPE, EKE
#histogram for energies
PE = .5*param['g']*param['rho']*h**2
PE = PE.reshape((tlen,param['ny'],param['nx']))[:,2:-2,2:-2].reshape((tlen,-1)).T
EPE = (PE.T - PE.mean(axis=1)).flatten()    # histogram from eddy potential energy
print('EPE done.')
del PE 

KE = .5*param['rho']*param['H']*(IuT.dot(u**2) + IvT.dot(v**2))
KE = KE.reshape((tlen,param['ny'],param['nx']))[:,2:-2,2:-2].reshape((tlen,-1)).T
EKE = (KE.T - KE.mean(axis=1)).flatten() # and from eddy kinetic energy
print('EKE done.')
del KE

#cut off outliers
EPEstd = EPE.std()
EKEstd = EKE.std()

s = 1.  # threshold
mask = np.logical_not((EPE > -s*EPEstd)*(EPE < 0.1*s*EPEstd)) + np.logical_not((EKE > -s*EKEstd)*(EKE < 0.1*s*EKEstd))
EPE = EPE[~mask]
EKE = EKE[~mask]
print('EPE EKE mask done.')

KEPEH,KE_edges,PE_edges = np.histogram2d(EKE,EPE,600)
print('EKE EPE histogram done.')
del EPE,EKE

KE_mid = KE_edges[:-1] + np.diff(KE_edges)/2.
PE_mid = PE_edges[:-1] + np.diff(PE_edges)/2.

## STORING
dic = dict()
all_var2export = ['KEPEH','KE_mid','PE_mid','KE_edges','PE_edges']

for v in all_var2export:
    exec('dic[v] ='+v)
    
np.save(runpath+'/analysis/KEPE_hist.npy',dic)

print('Everything stored.')