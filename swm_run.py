## SHALLOW WATER MODEL
from __future__ import print_function       # tested with python 3.6 and 2.7.12

# path
import os
path = os.getcwd() + '/'

# import modules
import numpy as np                          # version 1.11.3-py36
from scipy import sparse                    # version 0.19-py36
import time as tictoc
from netCDF4 import Dataset                 # version 1.2.4-py36, hdf5 version 1.8.17-py36, hdf4 version 4.2.12-py36
import glob
import zipfile

## import all functions
exec(open(path+'swm_param.py').read())
exec(open(path+'swm_operators.py').read())
exec(open(path+'swm_rhs.py').read())
exec(open(path+'swm_integration.py').read())
exec(open(path+'swm_output.py').read())

## set all parameters and initial conditions, and run model
u,v,eta = set_param()
u,v,eta = time_integration(u,v,eta)
