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
runfolder = [3,10]
print('Compare mean plots from run ' + str(runfolder))

## read data

runpath1 = path+'data/run%04i' % runfolder[0]
D1 = np.load(runpath1+'/analysis/mean.npy').all()
param1 = np.load(runpath1+'/param.npy').all()

runpath2 = path+'data/run%04i' % runfolder[1]
D2 = np.load(runpath2+'/analysis/mean.npy').all()
param2 = np.load(runpath2+'/param.npy').all()

runpath3 = path+'stoch/data/run%04i' % 13
D3 = np.load(runpath3+'/analysis/mean.npy').all()
param3 = np.load(runpath3+'/param.npy').all()

runpath4 = path+'stoch/data/run%04i' % 12
D4 = np.load(runpath4+'/analysis/mean.npy').all()
param4 = np.load(runpath4+'/param.npy').all()

runpath5 = path+'stoch/data/run%04i' % 14
D5 = np.load(runpath5+'/analysis/mean.npy').all()
param5 = np.load(runpath5+'/param.npy').all()

um,vm,hm = grid_interpolation(D2['um'],D2['vm'],D2['hm'],param2,param1)

##
# low resolution
stitle = r'Taylor Diagram: $\overline{u}$, $\overline{v}$, $\overline{h}$'

fig,ax = taylor_diag(hm, D1['hm'], norm=True,label=r'$h$',color='C0',marker='o')
taylor_diag(um, D1['um'], norm=True,label=r'$u$',add=True,ax=ax,color='C0',marker='o')
taylor_diag(vm, D1['vm'], norm=True,label=r'$v$',add=True,ax=ax,color='C0',marker='o')

# plus backscatter
taylor_diag(hm, D3['hm'], norm=True,label=r'$h$',add=True,ax=ax,color='C3')
taylor_diag(um, D3['um'], norm=True,label=r'$u$',add=True,ax=ax,color='C3')
taylor_diag(vm, D3['vm'], norm=True,label=r'$v$',add=True,ax=ax,color='C3')

# plus backscatter
taylor_diag(hm, D4['hm'], norm=True,label=r'$h$',add=True,ax=ax,color='C1')
taylor_diag(um, D4['um'], norm=True,label=r'$u$',add=True,ax=ax,color='C1')
taylor_diag(vm, D4['vm'], norm=True,label=r'$v$',add=True,ax=ax,color='C1')

# plus backscatter
taylor_diag(hm, D5['hm'], norm=True,label=r'$h$',add=True,ax=ax,color='C4')
taylor_diag(um, D5['um'], norm=True,label=r'$u$',add=True,ax=ax,color='C4')
taylor_diag(vm, D5['vm'], norm=True,label=r'$v$',add=True,ax=ax,color='C4',title=stitle)

ax.scatter(-np.ones(3),-np.ones(3),5,'C0',marker='o',label=r'Low resolution, $\Delta x = 30$km')
ax.scatter(-np.ones(3),-np.ones(3),5,'C2',marker='o',label=r'High resolution, $\Delta x = 7.5$km')
ax.scatter(-np.ones(3),-np.ones(3),5,'C3',marker='^',label=r'$\Delta x = 30$km + backscatter, $n_{diss} = %.2f$' % param3['n_diss'])
ax.scatter(-np.ones(3),-np.ones(3),5,'C1',marker='^',label=r'$\Delta x = 30$km + backscatter, $n_{diss} = %.2f$' % param4['n_diss'])
ax.scatter(-np.ones(3),-np.ones(3),5,'C4',marker='^',label=r'$\Delta x = 30$km + backscatter, $n_{diss} = %.2f$' % param5['n_diss'])

plt.legend(loc=1,scatterpoints=3)

plt.tight_layout()
plt.savefig(path+'compare/Taylor_uvh_mean.png',dpi=150)
plt.close(fig)