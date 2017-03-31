## PRODUCE MEAN CALCULATIONS AND EXPORT AS .NPY
from __future__ import print_function
path = '/home/mkloewer/python/swm/'
import os; os.chdir(path) # change working directory
import numpy as np
from scipy import sparse
import time as tictoc
from netCDF4 import Dataset
import glob

# OPTIONS
runfolder = [14]
print('Calculating means from run ' + str(runfolder))

## read data
for r,i in zip(runfolder,range(len(runfolder))):
    runpath = path+'stoch/data/run%04i' % r
    
    if i == 0:
        e = np.load(runpath+'/e_sub.npy')[10*365:,:,:]
        time = np.load(runpath+'/t_sub.npy')[10*365:]
        print('run %i read.' % r)

    else:
        e = np.concatenate((e,np.load(runpath+'/e_sub.npy')))
        time = np.hstack((time,np.load(runpath+'/t_sub.npy')))
        print('run %i read.' % r)

## create ouputfolder
try:
    os.mkdir(runpath+'/analysis')
except:
   pass
   
## U,V,H mean
em = e.mean(axis=0)
print('e mean done.')

## STORING
dic = dict()
all_var2export = ['em']

for v in all_var2export:
    exec('dic[v] ='+v)
    
np.save(runpath+'/analysis/mean_e.npy',dic)
print('Everything stored.')