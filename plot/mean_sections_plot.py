## PLOT STREAMFUNCTION
from __future__ import print_function
#path = '/home/mkloewer/python/swm/'
path = '/network/aopp/cirrus/pred/kloewer/swm_bf_cntrl/data/'
path2 = '/home/kloewer/git/swm/'
path3 = '/network/aopp/cirrus/pred/kloewer/swm_back_ronew/'
outpath = '/network/home/aopp/kloewer/swm/paperplot/'

import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt
from cmocean import cm

plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'

# OPTIONS
runfolder = [0,6,0,3,4]
print('Compare mean plots from run ' + str(runfolder))

## read data

runpath1 = path+'run%04i' % runfolder[0]
D1 = np.load(runpath1+'/analysis/mean.npy').all()
param1 = np.load(runpath1+'/param.npy').all()

runpath2 = path+'run%04i' % runfolder[1]
D2 = np.load(runpath2+'/analysis/mean.npy').all()
param2 = np.load(runpath2+'/param.npy').all()

runpath3 = path3+'run%04i' % runfolder[2]
D3 = np.load(runpath3+'/analysis/mean.npy').all()
param3 = np.load(runpath3+'/param.npy').all()

runpath4 = path3+'run%04i' % runfolder[3]
D4 = np.load(runpath4+'/analysis/mean.npy').all()
param4 = np.load(runpath4+'/param.npy').all()

runpath5 = path3+'run%04i' % runfolder[4]
D5 = np.load(runpath5+'/analysis/mean.npy').all()
param5 = np.load(runpath5+'/param.npy').all()

# functions
def h2mat(h,param):
    return h.reshape((param['ny'],param['nx']))

def u2mat(u,param):
    return u.reshape((param['ny'],param['nx']-1))

def v2mat(v,param):
    return v.reshape((param['ny']-1,param['nx']))

def q2mat(q,param):
    return q.reshape((param['ny']+1,param['nx']+1))

## Plotting
fig,axs = plt.subplots(1,3,sharex=True,figsize=(12,4))


xi = 10

axs[0].plot(param1['y_u']/1e3,u2mat(D1['um'],param1)[:,xi],"C0",lw=2)
axs[0].plot(param2['y_u']/1e3,u2mat(D2['um'],param2)[:,4*xi],"C2",lw=2)
#axs[0].plot(param3['y_u']/1e3,u2mat(D3['um'],param3)[:,xi],"C1",lw=1)
axs[0].plot(param4['y_u']/1e3,u2mat(D4['um'],param4)[:,xi],"C3",lw=1)
#axs[0].plot(param5['y_u']/1e3,u2mat(D5['um'],param5)[:,xi],"C5",lw=1)

axs[1].plot(param1['y_v']/1e3,v2mat(D1['vm'],param1)[:,xi],"C0",lw=2,label=r'Low resolution, $\Delta x = $30km')
axs[1].plot(param2['y_v']/1e3,v2mat(D2['vm'],param2)[:,4*xi],"C2",lw=2,label=r'High resolution, $\Delta x = $7.5km')
#axs[1].plot(param3['y_v']/1e3,v2mat(D3['vm'],param3)[:,xi],"C1",lw=1)
axs[1].plot(param4['y_v']/1e3,v2mat(D4['vm'],param4)[:,xi],"C3",lw=1,label=r'LR + moderate backscatter')
#axs[1].plot(param5['y_v']/1e3,v2mat(D5['vm'],param5)[:,xi],"C5",lw=1)

axs[2].plot(param1['y_T']/1e3,h2mat(D1['etam'],param1)[:,xi],"C0",lw=2)
axs[2].plot(param2['y_T']/1e3,h2mat(D2['etam'],param2)[:,4*xi],"C2",lw=2)
#axs[2].plot(param3['y_T']/1e3,h2mat(D3['etam'],param3)[:,xi],"C1",lw=1)
axs[2].plot(param4['y_T']/1e3,h2mat(D4['etam'],param4)[:,xi],"C3",lw=1)
#axs[2].plot(param5['y_T']/1e3,h2mat(D5['etam'],param5)[:,xi],"C5",lw=1)

axs[0].set_xlim(0,param1['Ly']/1e3)
axs[1].set_ylim(-0.4,.27)
axs[0].set_xlabel("y [km]"+" at x = {:d}km".format(int(param1['dx']/1e3*(xi+1))))
axs[1].set_xlabel("y [km]"+" at x = {:d}km".format(int(param1['dx']/1e3*(xi+.5))))
axs[2].set_xlabel("y [km]"+" at x = {:d}km".format(int(param1['dx']/1e3*(xi+.5))))

axs[0].set_ylabel("m/s")
axs[1].set_ylabel("m/s")
axs[2].set_ylabel("m")

axs[0].set_title('a',loc='left',fontweight='bold')
axs[1].set_title('b',loc='left',fontweight='bold')
axs[2].set_title('c',loc='left',fontweight='bold')

axs[0].set_title(r"zonal velocity $\bar{u}$")
axs[1].set_title(r"meridional velocity $\bar{v}$")
axs[2].set_title(r"surface elevation $\bar{\eta}$")

axs[1].legend(loc=2)

plt.tight_layout()
plt.savefig(outpath+'plots/mean_section.pdf')
plt.close(fig)
