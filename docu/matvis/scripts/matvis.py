## PLOTTING OPERATOR MATRICES
from __future__ import print_function
path = '/home/mkloewer/python/swm/'
import os; os.chdir(path) # change working directory
import numpy as np
from scipy import sparse
import time as tictoc
import matplotlib.pyplot as plt
from cmocean import cm

# import functions
exec(open(path+'swm_operators.py').read())
exec(open(path+'swm_output.py').read())
param = dict()
param['output'] = 0
param['dat_type'] = np.float32

param['dx'] = 1
param['dy'] = 1
param['nx'] = 3
param['ny'] = 3
param['NT'] = param['nx']*param['ny']
param['Nu'] = (param['nx']-1)*param['ny']
param['Nv'] = (param['ny']-1)*param['nx']
param['Nq'] = (param['nx']+1)*(param['ny']+1)

set_grad_mat()
set_interp_mat()
set_lapl_mat()
set_arakawa_mat()

##

#mnames = ['GTx','GTy','Gux','Guy','Gvx','Gvy','Gqy','Gqx']
#mnames = ['Lu','Lv','LT']
mnames = ['IuT','ITu','IvT','ITv','Iuv','Ivu','ITq','IqT','Iqu','Iuq','Iqv','Ivq']

for mname in mnames:
    exec('M = '+mname)
    levs = np.sort(np.array(list(set(M.data))))
    linlevs = np.arange(1,len(levs)+1)
    
    # replace data by linlevs
    idx = []
    for l in levs:
        idx.append(M.data == l)
    
    for i,r in zip(idx,range(len(idx))):
        M.data[i] = linlevs[r]
    
    M = M.todense()
    M = np.ma.masked_array(M,mask=(M == 0))
    
    aspectratio = M.shape[0]/M.shape[1]
    
    fig,ax = plt.subplots(1,1,figsize=(6,5*aspectratio))
    
    if len(levs) > 1:
        cmapd = cm.thermal.from_list('cmapd',plt.cm.jet(np.linspace(0,1,len(linlevs))),len(linlevs))
    else:
        cmapd = cm.thermal.from_list('cmapd',plt.cm.gray([0,1]),2)


    q = ax.matshow(M,cmap=cmapd,vmin=.5,vmax=linlevs.max()+.5)
    cbar = plt.colorbar(q,ax=ax,ticks=linlevs,drawedges=True)
    cbar.ax.set_yticklabels(levs)
    ax.set_xlabel(r'$\mathbf{'+mname[0]+'}^'+mname[2]+'_'+mname[1]+'$ for $n_x =$%i, $n_y =$%i' % (param['nx'],param['ny']),fontsize=20)
    plt.tight_layout()
    fig.savefig(path+'matvis/img/'+mname+'.png',dpi=150)
    plt.close(fig)
    #plt.show()