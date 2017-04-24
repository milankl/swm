## POWER INPUT EXIT COMPUTION AND PLOTTING
from __future__ import print_function

# path
import os
path = os.path.dirname(os.getcwd()) + '/'   # on level above
os.chdir(path)                              # change working directory

import numpy as np
from scipy import sparse

# OPTIONS
# several entries in the list concatenates the runs and stores the result in the last folder
runfolder = [0]     
print('Calculating power in/out from run ' + str(runfolder))

InPower_T = []
ExPower_T = []

## read data
for r in runfolder: # calculate each run separately
    runpath = path+'data/run%04i' % r
    
    u = np.load(runpath+'/u_sub.npy')
    v = np.load(runpath+'/v_sub.npy')
    eta = np.load(runpath+'/eta_sub.npy')
    t = np.load(runpath+'/t_sub.npy')
    print('run %i read.' % r)

    ## read param
    global param
    param = np.load(runpath+'/param.npy').all()
    
    # import functions
    exec(open(path+'swm_param.py').read())
    exec(open(path+'swm_operators.py').read())
    
    set_grad_mat()
    set_interp_mat()
    set_coriolis()
    set_forcing()
    
    tlen = len(t)
    ## create ouputfolder
    try:
        os.mkdir(runpath+'/analysis')
    except:
        pass
    
    ## reshape u,v
    u = u.reshape((tlen,param['Nu'])).T
    v = v.reshape((tlen,param['Nv'])).T
    h = eta.reshape((tlen,param['NT'])).T + param['H']
    print('Reshape done.')
    
    h_q = ITq.dot(h)
    h_u = ITu.dot(h)
    h_v = ITv.dot(h)
    print('h_u, h_v, h_q done.')
    
    ## input
    InPower = ((u.T*Fx/h_u.T).T*param['rho']*param['H']).mean(axis=1)
    print('Input Power done.')
    
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
    
    ExPower_u = (param['nu_B']*(u*bidiff_u)*param['rho']*param['H']).mean(axis=1)
    ExPower_v = (param['nu_B']*(v*bidiff_v)*param['rho']*param['H']).mean(axis=1)
    
    print('Exit Power done.')
    
    # Interpolation
    InPower_T.append(IuT.dot(InPower))
    ExPower_T.append(IuT.dot(ExPower_u) + IvT.dot(ExPower_v))

# Averaging over runs
InPower_T = np.array(InPower_T).mean(axis=0)
ExPower_T = np.array(ExPower_T).mean(axis=0)

## STORING
dic = dict()
all_var2export = ['InPower_T','ExPower_T']
for v in all_var2export:
    exec('dic[v] ='+v)
    
np.save(runpath+'/analysis/power_map.npy',dic)
print('Everything stored.')
