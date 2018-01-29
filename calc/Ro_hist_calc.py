## HISTOGRAM COMPUTATIONS FOR REYNOLDS AND ROSSBY NUMBERS
from __future__ import print_function

# path
import os
#path = '/network/aopp/cirrus/pred/kloewer/swm_bf_cntrl/'
path = '/network/aopp/cirrus/pred/kloewer/swm_back_ronew/'
os.chdir(path)                              # change working directory

# import functions
funpath = '/network/home/aopp/kloewer/git/swm/'
exec(open(funpath+'swm_param.py').read())
exec(open(funpath+'swm_operators.py').read())

import numpy as np
from scipy import sparse

# OPTIONS
runfolder = [5]
print('Calculating Rossby histograms from run ' + str(runfolder))

RoH = []
Ro_mean = []
Ro_median = []

## read data - calculate each run separately
for r,i in zip(runfolder,range(len(runfolder))):
    runpath = path+'run%04i' % r

    if i == 0:
        skip = 5*365
    else:
        skip = 0

    u = np.load(runpath+'/u_sub.npy')[skip:,...]
    v = np.load(runpath+'/v_sub.npy')[skip:,...]
    #eta = np.load(runpath+'/eta_sub.npy')[skip:,...]
    t = np.load(runpath+'/t_sub.npy')[skip:,...]
    print('run %i read.' % r)

    ## read param
    global param
    param = np.load(runpath+'/param.npy').all()

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

    ## COMPUTE ROSSBY VIA DEFORMATION RATE

    # Deformation rate / coriolis parameter
    Ro = (np.sqrt((Gux.dot(u) - Gvy.dot(v))**2 + IqT.dot((Guy.dot(u) + Gvx.dot(v))**2)).T / f_T).flatten()
    print('Ro computed.')

    Ro_mean.append(Ro.mean())
    Ro_median.append(np.median(Ro))
    Ro = np.log10(Ro)

    # histogram
    Ro_min = -5.    # in log scale
    Ro_max = 0.5
    N = 300

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

np.save(runpath+'/analysis/Ro_defo_hist.npy',dic)
