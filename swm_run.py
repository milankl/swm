## SHALLOW WATER MODEL
"""
Copyright (C) 2017,  Milan Kloewer (milankloewer@gmx.de)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

A documentation, in a preliminary version though, is available via
http://www.github.com/milankl/swm/docu/swm_documentation.pdf

"""

from __future__ import print_function       # tested with python 3.6 and 2.7.12

# path
import os
path = os.getcwd() + '/'         # this should be the path of this file

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
