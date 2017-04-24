## HISTOGRAM COMPUTATIONS FOR REYNOLDS AND ROSSBY NUMBERS
from __future__ import print_function

# path
import os
path = os.path.dirname(os.getcwd()) + '/'   # on level above
os.chdir(path)                              # change working directory

import numpy as np
from scipy import sparse

# OPTIONS
runfolder = [1]
print('Calculating Rossby histograms from run ' + str(runfolder))

RoH = []
Ro_mean = []
Ro_median = []

## read data - calculate each run separately
for r,i in zip(runfolder,range(len(runfolder))):
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
    
    tlen = len(t)
    ## create ouputfolder
    try:
        os.mkdir(runpath+'/analysis')
    except:
        pass
    
    ## reshape u,v
    u = u.reshape((tlen,param['Nu'])).T
    v = v.reshape((tlen,param['Nv'])).T
    print('Reshape done.')
        
    ## COMPUTE REYNOLDS, ROSSBY
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
    del u,v
    
    # actual number
    Ro = (adv_term / cor_term).flatten()
    print('Ro computed.')
    del adv_term, cor_term
    
    Ro_mean.append(Ro.mean())
    Ro_median.append(np.median(Ro))
    Ro = np.log10(Ro)
    
    # histogram
    Ro_min = -3     # in log scale
    Ro_max = 0
    N = 200
    
    RoH_temp,Ro_edges = np.histogram(Ro,np.linspace(Ro_min,Ro_max,N))
    print('Ro histogram done.')
    del Ro
    
    # store each run in a list
    RoH.append(RoH_temp)
    
    Ro_mid = Ro_edges[:-1] + np.diff(Ro_edges)/2.


RoH = np.array(RoH).sum(axis=0)
Ro_mean = np.array(Ro_mean).mean()
Ro_median = np.median(np.array(Ro_median))  #actually median of medians though...

## STORING in last 
dic = dict()
all_var2export = ['RoH','Ro_mid','Ro_edges','Ro_mean','Ro_median']

for vars in all_var2export:
    exec('dic[vars] ='+vars)
    
np.save(runpath+'/analysis/Ro_hist.npy',dic)