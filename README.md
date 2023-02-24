[![Project Status: Inactive – The project has reached a stable, usable state but is no longer being actively developed; support/maintenance will be provided as time allows.](https://www.repostatus.org/badges/latest/inactive.svg)](https://www.repostatus.org/#inactive)

# swm - Shallow water model
Shallow water equations solver with finite differences. Arakawa C-grid, Arakawa and Lamb advection scheme, Shchepetkin and O'Brian-like biharmonic diffusion operator. 4th order Runge-Kutta for time integration.

A documentation is available at http://www.github.com/milankl/swm/tree/master/docu

This model is written in python, relies on numpy, scipy and netCDF4.

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

# RUN

You may want to change settings in swm_param.py before you run the model. Then simply do

     python swm_run.py
     
# Energy budget-based backscatter

An implementation of the energy budget-based backscatter as described in Kloewer et al 2018, Ocean Modelling can be found in the separate branch "backscatter".

Copyright (C) 2018,  Milan Kloewer (milan.kloewer@physics.ox.ac.uk, milank.de)

# Citation

If you use this package in your research, please cite us

> Klöwer M, MF Jansen, M Claus, RJ Greatbatch and S Thomsen, 2018. *Energy budget-based backscatter in a shallow water model of a double gyre basin*, __Ocean Modelling__, 132, [10.1016/j.ocemod.2018.09.006](https://doi.org/10.1016/j.ocemod.2018.09.006)
