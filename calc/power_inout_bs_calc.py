## POWER INPUT EXIT CALCULATION
from __future__ import print_function

# path
import os
path = '/network/aopp/cirrus/pred/kloewer/swm_back_ronew/'
os.chdir(path)                              # change working directory

import numpy as np
from scipy import sparse

# OPTIONS
# several entries in the list concatenates the runs and stores the result in the last folder
runfolder = [4]
print('Calculating power in/out from run ' + str(runfolder))

InPower_T = []
ExPower_T = []
BfricPower_T = []
BackPower_T = []

## read data
for r in runfolder: # calculate each run separately
    runpath = path+'run%04i' % r

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
    funpath = '/network/home/aopp/kloewer/git/swm/'
    exec(open(funpath+'swm_param.py').read())
    exec(open(funpath+'swm_operators.py').read())

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
    e = e.reshape((tlen,param['NT'])).T
    print('Reshape done.')

    h_q = ITq.dot(h)
    h_u = ITu.dot(h)
    h_v = ITv.dot(h)
    print('h_u, h_v, h_q done.')

    ## input
    # Fx is actually Fx/rho
    InPower = param['rho']*(u.T*Fx).mean(axis=0)
    print('Input Power done.')

    # Shchepetkin and O'Brien divergence of a tensor formulation
    hS = ((Gux.dot(u)-Gvy.dot(v))*h,(G2vx.dot(v) + G2uy.dot(u))*h_q)
    diff_u = (GTx.dot(hS[0]) + Gqy.dot(hS[1])) / h_u
    diff_v = (Gqx.dot(hS[1]) - GTy.dot(hS[0])) / h_v

    # hS will be used for the backscatter term further down

    # biharmonic stress tensor R = (R11, R12, R12, -R11), store only R11, R12
    hR = ((Gux.dot(diff_u) - Gvy.dot(diff_v))*h, (G2vx.dot(diff_v) + G2uy.dot(diff_u))*h_q)


    del h_q, diff_u, diff_v

    # without the 1/h factor compared to the calculation in the rhs
    bidiff_u = (GTx.dot(hR[0]) + Gqy.dot(hR[1]))
    bidiff_v = (Gqx.dot(hR[1]) - GTy.dot(hR[0]))

    del hR

    ExPower_u = -param['nu_B']*param['rho']*(u*bidiff_u).mean(axis=1)
    ExPower_v = -param['nu_B']*param['rho']*(v*bidiff_v).mean(axis=1)

    del bidiff_u, bidiff_v

    print('Biharmonic dissipation term done.')

    BfricPower = -param['rho']*param['c_D']*((IuT.dot(u**2) + IvT.dot(v**2))**(3./2.)).mean(axis=1)

    print('Bottom friction dissipation term done.')

    # Backscatter term
    e_over_h = e/h
    nu_back = -param['c_back']*param['max_dxdy']*np.sqrt(2*e_over_h.clip(0,e_over_h.max()))
    del e_over_h, h, e

    nu_back_hS0 = nu_back*hS[0]
    nu_back_hS1 = ITq.dot(nu_back)*hS[1]

    del nu_back, hS

    # without the 1/h factor compared to the calculation in the rhs
    back_diff_u = (GTx.dot(nu_back_hS0) + Gqy.dot(nu_back_hS1))
    back_diff_v = (Gqx.dot(nu_back_hS1) - GTy.dot(nu_back_hS0))

    del nu_back_hS0, nu_back_hS1, h_u, h_v

    BackPower_u = param['rho']*(u*back_diff_u).mean(axis=1)
    BackPower_v = param['rho']*(v*back_diff_v).mean(axis=1)

    del back_diff_u, back_diff_v

    print('Backscatter term done.')

    # Interpolation
    InPower_T.append(IuT.dot(InPower))
    ExPower_T.append(IuT.dot(ExPower_u) + IvT.dot(ExPower_v))
    BfricPower_T.append(BfricPower)
    BackPower_T.append(IuT.dot(BackPower_u) + IvT.dot(BackPower_v))

# Averaging over runs
InPower_T = np.array(InPower_T).mean(axis=0)
ExPower_T = np.array(ExPower_T).mean(axis=0)
BfricPower_T = np.array(BfricPower_T).mean(axis=0)
BackPower_T = np.array(BackPower_T).mean(axis=0)

## STORING
dic = dict()
all_var2export = ['InPower_T','ExPower_T','BfricPower_T','BackPower_T']
for v in all_var2export:
    exec('dic[v] ='+v)

np.save(runpath+'/analysis/power_map.npy',dic)
print('Everything stored.')
