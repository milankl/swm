## TAYLOR DIAGRAM PLOT
from __future__ import print_function
path = '/home/mkloewer/python/swm/'
import os; os.chdir(path) # change working directory
import numpy as np
from scipy import sparse
import time as tictoc
from netCDF4 import Dataset
import glob
import matplotlib.pyplot as plt
from cmocean import cm

exec(open(path+'analyses/taylor_diagram.py').read())

# OPTIONS
runfolder = [3,10,11,12,13]
print('Compare mean plots from run ' + str(runfolder))

## read data

runpath1 = path+'data/run%04i' % runfolder[0]
D1 = np.load(runpath1+'/analysis/mean.npy').all()
param1 = np.load(runpath1+'/param.npy').all()

runpath2 = path+'data/run%04i' % runfolder[1]
D2 = np.load(runpath2+'/analysis/mean.npy').all()
param2 = np.load(runpath2+'/param.npy').all()

runpath3 = path+'data/run%04i' % runfolder[2]
D3 = np.load(runpath3+'/analysis/mean.npy').all()
param3 = np.load(runpath3+'/param.npy').all()

runpath4 = path+'data/run%04i' % runfolder[3]
D4 = np.load(runpath4+'/analysis/mean.npy').all()
param4 = np.load(runpath4+'/param.npy').all()

runpath5 = path+'data/run%04i' % runfolder[4]
D5 = np.load(runpath5+'/analysis/mean.npy').all()
param5 = np.load(runpath5+'/param.npy').all()

um128,vm128,hm128 = grid_interpolation(D2['um'],D2['vm'],D2['hm'],param2,param1)
um32,vm32,hm32 = grid_interpolation(D2['um'],D2['vm'],D2['hm'],param2,param3)
um64,vm64,hm64 = grid_interpolation(D2['um'],D2['vm'],D2['hm'],param2,param4)
um256,vm256,hm256 = grid_interpolation(D2['um'],D2['vm'],D2['hm'],param2,param5)

##
# low resolution
stitle = r'Taylor Diagram: $\overline{u}$, $\overline{v}$, $\overline{h}$'

fig,ax = taylor_diag(hm128, D1['hm'], norm=True,label=r'$h$',color='C0',marker='o')
taylor_diag(um128, D1['um'], norm=True,label=r'$u$',add=True,ax=ax,color='C0',marker='o')
taylor_diag(vm128, D1['vm'], norm=True,label=r'$v$',add=True,ax=ax,color='C0',marker='o')

# plus backscatter
taylor_diag(hm32, D3['hm'], norm=True,label=r'$h$',add=True,ax=ax,color='C3')
taylor_diag(um32, D3['um'], norm=True,label=r'$u$',add=True,ax=ax,color='C3')
taylor_diag(vm32, D3['vm'], norm=True,label=r'$v$',add=True,ax=ax,color='C3')

# plus backscatter
taylor_diag(hm64, D4['hm'], norm=True,label=r'$h$',add=True,ax=ax,color='C1')
taylor_diag(um64, D4['um'], norm=True,label=r'$u$',add=True,ax=ax,color='C1')
taylor_diag(vm64, D4['vm'], norm=True,label=r'$v$',add=True,ax=ax,color='C1')

# plus backscatter
taylor_diag(hm256, D5['hm'], norm=True,label=r'$h$',add=True,ax=ax,color='C4')
taylor_diag(um256, D5['um'], norm=True,label=r'$u$',add=True,ax=ax,color='C4')
taylor_diag(vm256, D5['vm'], norm=True,label=r'$v$',add=True,ax=ax,color='C4',title=stitle)

ax.scatter(-np.ones(3),-np.ones(3),5,'C3',marker='^',label=r'$\Delta x = 120$km')
ax.scatter(-np.ones(3),-np.ones(3),5,'C1',marker='^',label=r'$\Delta x = 60$km')
ax.scatter(-np.ones(3),-np.ones(3),5,'C0',marker='o',label=r'$\Delta x = 30$km')
ax.scatter(-np.ones(3),-np.ones(3),5,'C4',marker='^',label=r'$\Delta x = 15$km')
ax.scatter(-np.ones(3),-np.ones(3),5,'C2',marker='o',label=r'$\Delta x = 7.5$km')

plt.legend(loc=1,scatterpoints=3)

plt.tight_layout()
plt.savefig(path+'compare/Taylor_uvh_convergence.png',dpi=150)
plt.close(fig)