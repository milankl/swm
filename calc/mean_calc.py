## PRODUCE MEAN CALCULATIONS AND EXPORT AS .NPY
from __future__ import print_function
path = '/home/mkloewer/github/swm/'
import os; os.chdir(path) # change working directory
import numpy as np
from scipy import sparse
import time as tictoc
from netCDF4 import Dataset
import glob

# OPTIONS
runfolder = [2]
print('Calculating means from run ' + str(runfolder))

## read data
for r,i in zip(runfolder,range(len(runfolder))):
    runpath = path+'data/run%04i' % r
    
    if i == 0:
        u = np.load(runpath+'/u_sub.npy')[1*365:,:,:]
        v = np.load(runpath+'/v_sub.npy')[1*365:,:,:]
        h = np.load(runpath+'/h_sub.npy')[1*365:,:,:]
        time = np.load(runpath+'/t_sub.npy')[1*365:]
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


## U,V,H mean
um = u.mean(axis=1)  # temporal averages
vm = v.mean(axis=1)
hm = h.mean(axis=1)
print('u,v,h mean done.')

# MEAN KINETIC ENERGY
mke = .5*param['rho']*param['H']*(IuT.dot(um**2) + IvT.dot(vm**2))
# Mean Potential Energy
mpe = .5*param['g']*param['rho']*hm**2
print('MKE, MPE done.')

# VARIANCE
uvar = u.var(axis=1)  # temporal averages
vvar = v.var(axis=1)
epe = h.var(axis=1)*param['g']*param['rho']/2.   # eddy potential energy proportional to hvar
# Eddy kinetic energy
eke = .5*param['rho']*param['H']*(IuT.dot(uvar) + IvT.dot(vvar))
print('Variances, EPE, EKE done.')

# RELATIVE VORTICITY MEAN
zm = (Gvx.dot(v) - Guy.dot(u)).mean(axis=1)

# SPEED MEAN
speedm = np.sqrt(IuT.dot(u**2) + IvT.dot(v**2)).mean(axis=1)
print('Rel Vort, Speed done.')

# REYNOLDS, ROSSBY, EKMAN NUMBER MEAN
u_T = IuT.dot(u)
v_T = IvT.dot(v)
print('u,v interpolation done.')

#advective term
adv_u = u_T*Gux.dot(u) + v_T*IqT.dot(Guy.dot(u))
adv_v = u_T*IqT.dot(Gvx.dot(v)) + v_T*Gvy.dot(v)
del u_T,v_T
adv_term = np.sqrt(adv_u**2 + adv_v**2)
del adv_u, adv_v
adv_termm = adv_term.mean(axis=0)   # spatial average
adv_term = adv_term.mean(axis=1)    # temporal average
print('Advection term done.')

#coriolis term
cor_term = (f_T*np.sqrt(IuT.dot(u**2) + IvT.dot(v**2)).T).T
cor_termm = cor_term.mean(axis=0)   # spatial average
cor_term = cor_term.mean(axis=1)    # temporal average
print('Coriolis term done.')

#diffusive term
diff_u = param['nu_B']*LLu.dot(u)
diff_v = param['nu_B']*LLv.dot(v)
diff_term = np.sqrt(IuT.dot(diff_u**2) + IvT.dot(diff_v**2))
del diff_u,diff_v
diff_termm = diff_term.mean(axis=0)    # spatial average
diff_term = diff_term.mean(axis=1)     # temporal average
print('Diffusion term done.')

# temporal averages
Re = adv_term / diff_term
Ek = diff_term / cor_term
Ro = adv_term / cor_term

# spatial averages
Rem = adv_termm / diff_termm
Ekm = diff_termm / cor_termm
Rom = adv_termm / cor_termm
print('Nondim numbers done.')

## ENERGY AND ENSTROPHY TIME SERIES
PEm = (.5*param['g']*param['rho']*h**2).mean(axis=0)
print('Potential Energy done.')

KEm = (.5*param['rho']*param['H']*(IuT.dot(u**2) + IvT.dot(v**2))).mean(axis=0)
print('Kinetic Energy done.')

PV = (f_q + (Gvx.dot(v) - Guy.dot(u)).T).T / ITq.dot(h+param['H'])
Zm = (PV**2*ITq.dot(h+param['H'])).mean(axis=0)
del PV
print('Enstrophy done.')

## STORING
dic = dict()
all_var2export = ['um','vm','hm','mke','eke','mpe','epe']
all_var2export += ['uvar','vvar','zm','speedm']
all_var2export += ['Re','Ek','Ro','Rem','Ekm','Rom']
all_var2export += ['time','PEm','KEm','Zm']

for v in all_var2export:
    exec('dic[v] ='+v)
    
np.save(runpath+'/analysis/mean.npy',dic)
print('Everything stored.')