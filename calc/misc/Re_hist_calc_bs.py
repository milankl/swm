## HISTOGRAM COMPUTATIONS FOR REYNOLDS AND ROSSBY NUMBERS
from __future__ import print_function

# path
import os
path = os.path.dirname(os.getcwd()) + '/'   # on level above
os.chdir(path)                              # change working directory

import numpy as np
from scipy import sparse

# OPTIONS
runfolder = [10]
print('Calculating Reynolds histograms from run ' + str(runfolder))

ReH = []
Re_mean = []
Re_median = []

## read data - calculate each run separately
for r,i in zip(runfolder,range(len(runfolder))):
    runpath = path+'data/run%04i' % r
    
    skip = 5*365
    u = np.load(runpath+'/u_sub.npy')[skip:,...]
    v = np.load(runpath+'/v_sub.npy')[skip:,...]
    eta = np.load(runpath+'/eta_sub.npy')[skip:,...]
    e = np.load(runpath+'/e_sub.npy')[skip:,...]
    t = np.load(runpath+'/t_sub.npy')[skip:,...]
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
    h = eta.reshape((tlen,param['NT'])).T + param['H']
    e = e.reshape((tlen,param['NT'])).T
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
    
    #diffusive term
    S = (Gux.dot(u)-Gvy.dot(v),G2vx.dot(v) + G2uy.dot(u))
    del u,v
    hS = (h*S[0],ITq.dot(h)*S[1])
    del S
    print('Stress tensor S done.')

    diff_u = (GTx.dot(hS[0]) + Gqy.dot(hS[1])) / ITu.dot(h)
    diff_v = (Gqx.dot(hS[1]) - GTy.dot(hS[0])) / ITv.dot(h)
    print('Harmonic part done.')

    # biharmonic stress tensor R = (R11, R12, R12, -R11), store only R11, R12
    R = (Gux.dot(diff_u) - Gvy.dot(diff_v), G2vx.dot(diff_v) + G2uy.dot(diff_u))
    del diff_u, diff_v
    nuhR = (param['nu_B']*h*R[0],param['nu_B']*ITq.dot(h)*R[1])
    del R
    print('Stress tensor R done.')
    
    bidiff_u = (GTx.dot(nuhR[0]) + Gqy.dot(nuhR[1])) / ITu.dot(h)
    bidiff_v = (Gqx.dot(nuhR[1]) - GTy.dot(nuhR[0])) / ITv.dot(h)

    del nuhR
    print('Biharmonic part done.')
    
    # backscatter
    e_over_h = e/h
    nu_back = -param['c_back']*param['max_dxdy']*np.sqrt(2*e_over_h.clip(0,e_over_h.max()))
    del e_over_h,e
    
    nu_back_hS0 = nu_back*hS[0]
    nu_back_hS1 = ITq.dot(nu_back)*hS[1]
    print('nu_back calculated.')
    del nu_back
    
    back_diff_u = (GTx.dot(nu_back_hS0) + Gqy.dot(nu_back_hS1)) / ITu.dot(h)
    back_diff_v = (Gqx.dot(nu_back_hS1) - GTy.dot(nu_back_hS0)) / ITv.dot(h)
    del nu_back_hS0,nu_back_hS1
    
    diff_term = np.sqrt(IuT.dot((-bidiff_u + back_diff_u)**2) + IvT.dot((-bidiff_v + back_diff_v)**2))
    #diff_term = np.sqrt(IuT.dot(bidiff_u**2) + IvT.dot(bidiff_v**2))
    print('Diff term done.')
    del bidiff_u,bidiff_v,back_diff_u,back_diff_v
    
    # actual number
    Re = (adv_term / diff_term).flatten()
    print('Re computed.')
    del adv_term, diff_term
    
    Re_mean.append(Re.mean())
    Re_median.append(np.median(Re))
    Re = np.log10(Re)
    
    # histogram
    Re_min = -3.    # in log scale
    Re_max = 5.
    N = 300
    
    ReH_temp,Re_edges = np.histogram(Re,np.linspace(Re_min,Re_max,N))
    print('Re histogram done.')
    del Re
    
    # store each run in a list
    ReH.append(ReH_temp)
    
    Re_mid = Re_edges[:-1] + np.diff(Re_edges)/2.


ReH = np.array(ReH).sum(axis=0)
Re_mean = np.array(Re_mean).mean()
Re_median = np.median(np.array(Re_median))  #actually median of medians though...

## STORING in last 
dic = dict()
all_var2export = ['ReH','Re_mid','Re_edges','Re_mean','Re_median']

for vars in all_var2export:
    exec('dic[vars] ='+vars)
    
np.save(runpath+'/analysis/Re_hist.npy',dic)
