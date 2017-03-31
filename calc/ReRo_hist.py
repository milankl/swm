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
runfolder = [12]
print('Calculating histograms from run ' + str(runfolder))

## read data
for r,i in zip(runfolder,range(len(runfolder))):
    runpath = path+'stoch/data/run%04i' % r
    
    if i == 0:
        u = np.load(runpath+'/u_sub.npy')[10*365:,:,:]
        v = np.load(runpath+'/v_sub.npy')[10*365:,:,:]
        h = np.load(runpath+'/h_sub.npy')[10*365:,:,:]        
        time = np.load(runpath+'/t_sub.npy')[10*365:]
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

#diffusive term
diff_u = param['nu_B']*LLu.dot(u)
diff_v = param['nu_B']*LLv.dot(v)
diff_term = np.sqrt(IuT.dot(diff_u**2) + IvT.dot(diff_v**2))
del diff_u,diff_v
print('Diffusion term done.')

# actual numbers
Re = np.log10(adv_term / diff_term).flatten()
Ro = np.log10(adv_term / cor_term).flatten()
del adv_term, cor_term, diff_term
print('Non dim numbers computed.')

# histogram
ReRoH,Re_edges,Ro_edges = np.histogram2d(Re,Ro,200)
print('ReRo histogram done.')
del Re,Ro

Re_mid = Re_edges[:-1] + np.diff(Re_edges)/2.
Ro_mid = Ro_edges[:-1] + np.diff(Ro_edges)/2.

## STORING
dic = dict()
all_var2export = ['ReRoH','Re_mid','Ro_mid','Re_edges','Ro_edges']

for v in all_var2export:
    exec('dic[v] ='+v)
    
np.save(runpath+'/analysis/ReRo_hist.npy',dic)