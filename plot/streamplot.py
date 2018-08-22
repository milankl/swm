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
runfolder = [0,6,3]
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
print(param3['n_diss'])
    
## CALCULATE
u = [0]*3
v = [0]*3
speed = [0]*3
xx = [0]*3
yy = [0]*3


global param
param = np.load(runpath1+'/param.npy').all()
exec(open(path2+'swm_operators.py').read())
exec(open(path2+'swm_output.py').read())
param['output'] = 0
set_grad_mat()
set_interp_mat()

u[0] = h2mat(IuT.dot(D1['um']))
v[0] = h2mat(IvT.dot(D1['vm']))
speed[0] = np.sqrt(u[0]**2 + v[0]**2)
xx[0],yy[0] = np.meshgrid(param['x_T'],param['y_T'])

# next run
param = np.load(runpath2+'/param.npy').all()
exec(open(path2+'swm_operators.py').read())
exec(open(path2+'swm_output.py').read())
param['output'] = 0
set_grad_mat()
set_interp_mat()

u[1] = h2mat(IuT.dot(D2['um']))
v[1] = h2mat(IvT.dot(D2['vm']))
speed[1] = np.sqrt(u[1]**2 + v[1]**2)

xx[1],yy[1] = np.meshgrid(param['x_T'],param['y_T'])

# next run
param = np.load(runpath3+'/param.npy').all()
exec(open(path2+'swm_operators.py').read())
exec(open(path2+'swm_output.py').read())
param['output'] = 0
set_grad_mat()
set_interp_mat()

u[2] = h2mat(IuT.dot(D3['um']))
v[2] = h2mat(IvT.dot(D3['vm']))
speed[2] = np.sqrt(u[2]**2 + v[2]**2)

xx[2],yy[2] = np.meshgrid(param['x_T'],param['y_T'])

# line width
lww = [4*np.sqrt(speed[0])/np.sqrt(speed[0].max())]
lww.append(4*np.sqrt(speed[1])/np.sqrt(speed[0].max()))
lww.append(4*np.sqrt(speed[2])/np.sqrt(speed[0].max()))

# wind profile
Fx = param['Fx0']*(np.cos(2*np.pi*(param['y_u']-param['Ly']/2)/param['Ly']) + 2*np.sin(np.pi*(param['y_u'] - param['Ly']/2)/param['Ly']))

## Plotting
fig,axs = plt.subplots(1,3,sharex=True,sharey=True,figsize=(12,5))
fig.tight_layout(rect=[0.02,.15,0.9,0.92])
fig.subplots_adjust(wspace=0.03,hspace=0.03)

pos = axs[0].get_position()
pos2 = axs[2].get_position()
cax = fig.add_axes([pos.x0,0.1,pos2.x1-pos.x0,0.03])

s = 1e3 # scaling factor m -> km

# wind profile
wax = fig.add_axes([0.888,pos.y0,0.1,pos.y1-pos.y0])
wax.set_ylim(0,param['Ly']/s)
wax.set_yticks([0,1000,2000,3000])
wax.set_yticklabels([])
wax.plot(Fx,param['y_u']/s)
wax.plot([0,0],[0,param['Ly']/s],'grey',lw=0.5)
#wax.set_title('d',fontweight='bold',loc='left')
wax.set_title(r'Wind stress $\tau$',loc='left')


wax.set_xlabel('[Pa]')
wax.set_xticks([-0.3,0,0.2])
wax.set_xticklabels([-0.3,0,0.2])

#plt.text(0, 1.1,r'Wind stress $\tau$',fontsize=13,ha='left',transform=wax.transAxes)

for j,i in zip([0,2,1],range(3)):
    strm = axs[j].streamplot(xx[i]/s,yy[i]/s,u[i],v[i],color=speed[i],density=2.5,linewidth=lww[i],arrowstyle='->')

cbar = fig.colorbar(strm.lines,cax=cax,orientation='horizontal')
cbar.set_label('Speed [m/s]')
axs[0].set_xticks([0,1000,2000,3000])
axs[0].set_yticks([0,1000,2000,3000])

axs[0].set_xlim(0,param['Lx']/s)
axs[0].set_ylim(0,param['Ly']/s)
axs[0].set_ylabel(r'$y$ [km]')
axs[1].set_xlabel(r'$x$ [km]')

axs[0].set_title(r'Low resolution, $\Delta x = $30km')
axs[1].set_title(r'LR + moderate backscatter')
axs[2].set_title(r'High resolution, $\Delta x = $7.5km')

axs[0].set_title('a',fontweight='bold',loc='left')
axs[1].set_title('b',fontweight='bold',loc='left')
axs[2].set_title('c',fontweight='bold',loc='left')

#plt.suptitle('Streamlines of mean circulation ($\overline{u},\overline{v}$)',x=0.186,y=0.98)

#plt.text(0,1.1,'Streamlines of mean circulation ($\overline{u},\overline{v}$)',fontsize=13,ha='left',transform=axs[0].transAxes)

#plt.suptitle(r'wind stress $\tau$',x=0.86,y=0.98)
plt.savefig(outpath+'plots/streamplot.eps')
plt.close(fig)
