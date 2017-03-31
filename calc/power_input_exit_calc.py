## POWER INPUT EXIT COMPUTION AND PLOTTING
from __future__ import print_function
path = '/home/mkloewer/python/swm/'
import os; os.chdir(path) # change working directory
import numpy as np
from scipy import sparse
import time as tictoc
from netCDF4 import Dataset
import glob

from scipy.integrate import cumtrapz

# OPTIONS
runfolder = [2,3]
print('Calculating power input exit from run ' + str(runfolder))

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

try:
    param['nu_B'] = param['B']
except:
    pass

set_grad_mat()
set_interp_mat()
set_lapl_mat()
set_coriolis()
set_forcing()

tlen = len(time)
## create ouputfolder
try:
    os.mkdir(runpath+'/analysis')
except:
    pass
   
## reshape u,v
u = u.reshape((tlen,param['Nu'])).T
v = v.reshape((tlen,param['Nv'])).T
h = h.reshape((tlen,param['NT'])).T + param['H']
print('Reshape done.')

h_q = ITq.dot(h)
h_u = ITu.dot(h)
h_v = ITv.dot(h)
print('h_u, h_v, h_q done.')

## input
InPower = ((u.T*Fx).T*param['rho']*param['H']).mean(axis=0)
InEnergy = cumtrapz(InPower,time,initial=0)

print('Input Power and Energy done.')

# Shchepetkin and O'Brien divergence of a tensor formulation
hS = ((Gux.dot(u)-Gvy.dot(v))*h,(G2vx.dot(v) + G2uy.dot(u))*h_q)
diff_u = (GTx.dot(hS[0]) + Gqy.dot(hS[1])) / h_u
diff_v = (Gqx.dot(hS[1]) - GTy.dot(hS[0])) / h_v

del hS

# biharmonic stress tensor R = (R11, R12, R12, -R11), store only R11, R12
hR = ((Gux.dot(diff_u) - Gvy.dot(diff_v))*h, (G2vx.dot(diff_v) + G2uy.dot(diff_u))*h_q)

del h_q, diff_u, diff_v

bidiff_u = (GTx.dot(hR[0]) + Gqy.dot(hR[1])) / h_u
bidiff_v = (Gqx.dot(hR[1]) - GTy.dot(hR[0])) / h_v

del hR, h_u, h_v

print('Biharmonic dissipation term done.')

ExPower_u = (param['nu_B']*(u*bidiff_u)*param['rho']*param['H']).mean(axis=0)
ExPower_v = (param['nu_B']*(v*bidiff_v)*param['rho']*param['H']).mean(axis=0)

ExEnergy = cumtrapz(ExPower_u,time,initial=0) + cumtrapz(ExPower_v,time,initial=0)

print('Exit Energy and Power done.')
## STORING
dic = dict()
all_var2export = ['time','InPower','InEnergy','ExPower_u','ExPower_v','ExEnergy']
for v in all_var2export:
    exec('dic[v] ='+v)
    
np.save(runpath+'/analysis/power_energy_timeseries.npy',dic)
print('Everything stored.')
