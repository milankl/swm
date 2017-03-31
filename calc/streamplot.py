## PLOT STREAMFUNCTION
from __future__ import print_function
path = '/home/mkloewer/python/swm/'
import os; os.chdir(path) # change working directory
import numpy as np
from scipy import sparse
from scipy.integrate import cumtrapz
import matplotlib.pyplot as plt
import time as tictoc
from netCDF4 import Dataset
import glob
from cmocean import cm

# OPTIONS
runfolder = [3,10]
print('Compare mean plots from run ' + str(runfolder))

## read data

runpath1 = path+'data/run%04i' % runfolder[0]
D1 = np.load(runpath1+'/analysis/mean.npy').all()
param1 = np.load(runpath1+'/param.npy').all()

runpath2 = path+'data/run%04i' % runfolder[1]
D2 = np.load(runpath2+'/analysis/mean.npy').all()
param2 = np.load(runpath2+'/param.npy').all()

runpath3 = path+'stoch/data/run%04i' % 12
D3 = np.load(runpath3+'/analysis/mean.npy').all()
param3 = np.load(runpath3+'/param.npy').all()
    
## CALCULATE
u = [0]*3
v = [0]*3
speed = [0]*3
xx = [0]*3
yy = [0]*3


global param
param = np.load(runpath1+'/param.npy').all()
exec(open(path+'swm_operators.py').read())
exec(open(path+'swm_output.py').read())
param['output'] = 0
set_grad_mat()
set_interp_mat()

u[0] = h2mat(IuT.dot(D1['um']))
v[0] = h2mat(IvT.dot(D1['vm']))
speed[0] = np.sqrt(u[0]**2 + v[0]**2)
xx[0],yy[0] = np.meshgrid(param['x_T'],param['y_T'])

# next run
param = np.load(runpath2+'/param.npy').all()
exec(open(path+'swm_operators.py').read())
exec(open(path+'swm_output.py').read())
param['output'] = 0
set_grad_mat()
set_interp_mat()

u[1] = h2mat(IuT.dot(D2['um']))
v[1] = h2mat(IvT.dot(D2['vm']))
speed[1] = np.sqrt(u[1]**2 + v[1]**2)

xx[1],yy[1] = np.meshgrid(param['x_T'],param['y_T'])

# next run
param = np.load(runpath3+'/param.npy').all()
exec(open(path+'swm_operators.py').read())
exec(open(path+'swm_output.py').read())
param['output'] = 0
set_grad_mat()
set_interp_mat()

u[2] = h2mat(IuT.dot(D3['um']))
v[2] = h2mat(IvT.dot(D3['vm']))
speed[2] = np.sqrt(u[2]**2 + v[2]**2)

xx[2],yy[2] = np.meshgrid(param['x_T'],param['y_T'])

# line width
lww = [3*np.sqrt(speed[0])/np.sqrt(speed[0].max())]
lww.append(3*np.sqrt(speed[1])/np.sqrt(speed[0].max()))
lww.append(3*np.sqrt(speed[2])/np.sqrt(speed[0].max()))


## Plotting
fig,axs = plt.subplots(1,3,sharex=True,sharey=True,figsize=(12,6))
fig.tight_layout(rect=[0,.1,1,0.95])
fig.subplots_adjust(wspace=0.03,hspace=0.03)

pos = axs[0].get_position()
pos2 = axs[2].get_position()
cax = fig.add_axes([pos.x0,0.1,pos2.x1-pos.x0,0.03])

for j,i in zip([0,2,1],range(3)):
    strm = axs[j].streamplot(xx[i],yy[i],u[i],v[i],color=speed[i],density=2,linewidth=lww[i])

cbar = fig.colorbar(strm.lines,cax=cax,orientation='horizontal')
cbar.set_label('Speed [m/s]')
axs[0].set_xticks([])
axs[0].set_yticks([])

axs[0].set_xlim(0,param['Lx'])
axs[0].set_ylim(0,param['Ly'])

axs[0].set_title(r'Low resolution, $\Delta x = 30$km')
axs[1].set_title(r'Low resolution + backscatter')
axs[2].set_title(r'High resolution, $\Delta x = 7.5$km')


plt.suptitle('Streamlines of ($\overline{u},\overline{v}$)')
plt.savefig(path+'compare/streamplot.png',dpi=150)
plt.close(fig)