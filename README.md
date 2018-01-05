# Shallow water model
Shallow water equations solver with finite differences. Arakawa C-grid, Arakawa and Lamb advection scheme, Shchepetkin and O'Brian-like biharmonic diffusion operator. 4th order Runge-Kutta for time integration.

A documentation is available at http://www.github.com/milankl/swm/docu

This model is written in python, relies on numpy, scipy and netCDF4 as well as parallel_sparsetools (see folder) for a parallel matrix-vector multiplication.

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

# HOW TO INSTALL

For the parallel matrix-vector multiplication from parallel_sparsetools, go to that folder and do:
     
     python setup.py install
     
Otherwise comment the line 

     import _parallel_sparsetools
     
in swm_run.py, and

     os.environ['OMP_NUM_THREADS'] = str(1)
     sparse.csr_matrix._mul_vector = _mul_vector

in swm_param.py. parallel_sparsetools might yield a speed-up of up to 2.5x for up to 4 cores on some machines.

# RUN

You may want to change settings in swm_param.py before you run the model. Then simply do

     python swm_run.py

Copyright (C) 2017,  Milan Kloewer (milan.kloewer@physics.ox.ac.uk, milank.de)
