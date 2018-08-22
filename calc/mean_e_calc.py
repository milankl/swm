## PRODUCE MEAN CALCULATIONS AND EXPORT AS .NPY
from __future__ import print_function
path = '/home/mkloewer/python/swm/'
import os; os.chdir(path) # change working directory
import numpy as np
from scipy import sparse
import time as tictoc
from netCDF4 import Dataset

# OPTIONS
runfolder = 15
print('Calculating subgrid-EKE means from run ' + str(runfolder))

## read data
runpath = path+'data/run%04i' % runfolder

skip = 5*365
e = np.load(runpath+'/e_sub.npy')[skip:,:,:]
print('run %i read.' % runfolder)

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
