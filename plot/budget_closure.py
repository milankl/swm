from __future__ import print_function

path = '/network/aopp/cirrus/pred/kloewer/swm_bf_cntrl/data/'
dpath = '/network/aopp/cirrus/pred/kloewer/swm_back_ronew/'
outpath = '/network/home/aopp/kloewer/swm/paperplot/'

import os; os.chdir(path) # change working directory
import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt
import time as tictoc
from netCDF4 import Dataset
import glob
from matplotlib.colors import BoundaryNorm,LogNorm
import cmocean

plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'

# functions
def h2mat(h,param):
    return h.reshape((param['ny'],param['nx']))

def u2mat(u,param):
    return u.reshape((param['ny'],param['nx']-1))

def v2mat(v,param):
    return v.reshape((param['ny']-1,param['nx']))

def q2mat(q,param):
    return q.reshape((param['ny']+1,param['nx']+1))

# OPTIONS
runfolders_cntrl = [12,10,0,11,6,9]
runfolders_back = [0,1,2,3,6,7,8,4]

runfolders_start = np.array([1,0,5,2,5,0])*365               # start year for computations of energy budget

budget_terms_cntrl = np.zeros((len(runfolders_cntrl),4))     # four terms for cntrl: dE/dt, wind, bottom friction, viscosity
budget_terms_back = np.zeros((len(runfolders_back),5))       # five terms for back:  dE/dt, wind, bottom friction, viscosity, backscatter

res_cntrl = np.zeros(len(runfolders_cntrl))
dx_cntrl = np.zeros(len(runfolders_cntrl))
rdiss = np.zeros(len(runfolders_back))

# CONTROL
for i,r in enumerate(runfolders_cntrl):

    runpath = path+'run%04i' % r
    D = np.load(runpath+'/analysis/power_map.npy').all()
    E = np.load(runpath+'/analysis/mean_timeseries.npy').all()
    param = np.load(runpath+'/param.npy').all()

    Estart = (E['KEm']+E['PEm'])[runfolders_start[i]]
    Eend = (E['KEm']+E['PEm'])[-1]
    dt = 4*param['output_dt']*(len(E['KEm'][runfolders_start[i]:])-1)

    budget_terms_cntrl[i,0] = (Eend - Estart)/dt  
    budget_terms_cntrl[i,1] = D['InPower_T'].mean()
    budget_terms_cntrl[i,2] = D['BfricPower_T'].mean()
    budget_terms_cntrl[i,3] = D['ExPower_T'].mean()
    
    res_cntrl[i] = param['nx']
    dx_cntrl[i] = param['dx']
    
# BACKSCATTER    
for i,r in enumerate(runfolders_back):

    runpath = dpath+'run%04i' % r       # use different path
    D = np.load(runpath+'/analysis/power_map.npy').all()
    E = np.load(runpath+'/analysis/mean_timeseries.npy').all()
    param = np.load(runpath+'/param.npy').all()
    
    Estart = (E['KEm']+E['PEm'])[5*365]     # for backscatter runs discard first 5 years
    Eend = (E['KEm']+E['PEm'])[-1]
    dt = 4*param['output_dt']*(len(E['KEm'][5*365:])-1)
    
    budget_terms_back[i,0] = (Eend - Estart)/dt
    budget_terms_back[i,1] = D['InPower_T'].mean()
    budget_terms_back[i,2] = D['BfricPower_T'].mean()
    budget_terms_back[i,3] = D['ExPower_T'].mean()
    budget_terms_back[i,4] = D['BackPower_T'].mean()

    rdiss[i] = param['n_diss']


closure_cntrl = budget_terms_cntrl.sum(axis=1)
closure_back = budget_terms_back.sum(axis=1)

closure_norm_cntrl = np.sqrt((budget_terms_cntrl**2).sum(axis=1))
closure_norm_back = np.sqrt((budget_terms_back**2).sum(axis=1))

# treat viscosity and backscatter as one
budget_terms_back[:,3] = budget_terms_back[:,3:].sum(axis=1)
budget_terms_back[:,4] = 0

dissipation = budget_terms_cntrl[:,2:].sum(axis=1)

## PLOTTING 1 versus backscatter strength Rdiss

s = 1e3
p = 0.001

norm = (closure_norm_back.mean() + closure_norm_cntrl.mean())/2*s*p

fig,(ax1,ax) = plt.subplots(2,1,sharex=True)

for i in range(len(runfolders_cntrl)):
    ax1.plot(0,budget_terms_cntrl[i,0]*s,"C"+str(i)+"x",alpha=.6)
    ax1.plot(0,budget_terms_cntrl[i,1]*s,"C"+str(i)+"s",alpha=.6)
    ax1.plot(0,budget_terms_cntrl[i,2]*s,"C"+str(i)+"^",alpha=.6)
    ax1.plot(0,budget_terms_cntrl[i,3]*s,"C"+str(i)+"o",alpha=.6)

ax1.plot(-10,0,"x",color="grey",label="dE/dt")
ax1.plot(-10,0,"s",color="grey",label="wind stress")
ax1.plot(-10,0,"^",color="grey",label="bottom friction")
ax1.plot(-10,0,"o",color="grey",label="viscosity + backscatter")


ax1.legend(loc=3,ncol=4,fontsize=6)

for i in range(len(runfolders_back)):
    ax1.plot(i+1,budget_terms_back[i,0]*s,"kx",alpha=.7)
    ax1.plot(i+1,budget_terms_back[i,1]*s,"ks",alpha=.7)
    ax1.plot(i+1,budget_terms_back[i,2]*s,"k^",alpha=.7)
    ax1.plot(i+1,budget_terms_back[i,3]*s,"ko",alpha=.7)

ax.plot([-1,10],[0,0],"grey",lw=0.2)
ax1.plot([-1,10],[0,0],"grey",lw=0.2)



lines = [0]*len(closure_cntrl)

for i,cc in enumerate(closure_cntrl):
    lines[i], = ax.plot(0,cc*s,"C"+str(i)+"+",ms=3,label=r"$\Delta x$ = {:.2f}km".format(dx_cntrl[i]/1e3))

first_legend = plt.legend(loc=3,handles=lines,title="Control runs",fontsize=6)
fl = plt.gca().add_artist(first_legend)

linef = ax.fill_between([-1,10],[-norm,-norm],[norm,norm],alpha=.2,label="budget unclosed by <0.1%")
line, = ax.plot(1+np.arange(len(rdiss)),closure_back*s,"k*",label=r"$\Delta x$ = 30km",alpha=.7)
plt.legend(loc=8,handles=[line,linef],title="Backscatter runs",fontsize=6)

ax1.set_ylim(-10,10)
ax.set_xlim(-1,len(rdiss)+1)
ax.set_xticks(np.arange(len(rdiss)+1))
ax.set_xticklabels([0,1,2,6,8,16,32,64,r"$\infty$"])

ax.set_xlabel(r"Backscatter strength $R_{diss}$")
ax.set_ylabel(r"[mW m$^{-2}$]")
ax1.set_ylabel(r"[mW m$^{-2}$]")
ax1.set_title("Energy budget")
ax.set_title("Budget closure")

ax1.set_title("a",loc="left",fontweight="bold")
ax.set_title("b",loc="left",fontweight="bold")


plt.tight_layout()
plt.savefig(outpath+'plots/budget_closure.pdf')
plt.close(fig)

## PLOTTING 2 - versus RESOLUTION

fig,(ax,ax1) = plt.subplots(2,1)

for i in range(len(runfolders_cntrl)):
    ax.plot(res_cntrl[i],budget_terms_cntrl[i,0]*s,"C"+str(i)+"x",alpha=.6)
    ax.plot(res_cntrl[i],budget_terms_cntrl[i,1]*s,"C"+str(i)+"s",alpha=.6)
    ax.plot(res_cntrl[i],budget_terms_cntrl[i,2]*s,"C"+str(i)+"^",alpha=.6)
    ax.plot(res_cntrl[i],budget_terms_cntrl[i,3]*s,"C"+str(i)+"o",alpha=.6)

ax.plot(-10,0,"x",color="grey",label="dE/dt")
ax.plot(-10,0,"s",color="grey",label="wind stress")
ax.plot(-10,0,"^",color="grey",label="bottom friction")
ax.plot(-10,0,"o",color="grey",label="viscosity")

for i in range(len(runfolders_cntrl)):
    ax1.plot(res_cntrl[i],budget_terms_cntrl[i,2]/dissipation[i]*100,"C"+str(i)+"^",alpha=.6)
    ax1.plot(res_cntrl[i],budget_terms_cntrl[i,3]/dissipation[i]*100,"C"+str(i)+"o",alpha=.6)


ax.plot(res_cntrl,budget_terms_cntrl[:,1]*s,"k",zorder=1,lw=.8)
ax.plot(res_cntrl,budget_terms_cntrl[:,2]*s,"k",zorder=1,lw=.8)
ax.plot(res_cntrl,budget_terms_cntrl[:,3]*s,"k",zorder=1,lw=.8)

ax1.plot(res_cntrl,budget_terms_cntrl[:,2]/dissipation*100,"k",zorder=1,lw=1)
ax1.plot(res_cntrl,budget_terms_cntrl[:,3]/dissipation*100,"k",zorder=1,lw=1)


ax.legend(loc=4,ncol=4,fontsize=6)
ax.plot([-1,1300],[0,0],"--",lw=1,color="grey")
ax1.plot([-1,1300],[0,0],"--",lw=1,color="grey")
ax1.plot([-1,1300],[100,100],"--",lw=1,color="grey")

ax1.plot(-100,0,"^",color="grey",label="bottom friction")
ax1.plot(-100,0,"o",color="grey",label="viscosity")
ax1.legend(loc=7)

ax.set_ylim(-10,10)
ax1.set_ylim(-5,105)
ax.set_xlim(-1,1050)
ax1.set_xlim(-1,1050)
ax1.set_xticks(res_cntrl)
ax.set_xticks(res_cntrl)
ax.set_xticklabels([120,60,30,15,7.5,3.75],rotation=90,fontsize=8)
ax1.set_xticklabels([r"32$^2$",r"64$^2$",r"128$^2$",r"256$^2$",r"512$^2$",r"1024$^2$"],fontsize=8)

ax.set_xlabel(r"Resolution $\Delta x$ [km]",fontsize=8)
ax1.set_xlabel(r"Grid cells $N^2$",fontsize=8)
ax.set_ylabel(r"[mW m$^{-2}$]")
ax1.set_ylabel("[%]")
ax.set_title("Energy budget")
ax1.set_title("Dissipation: Bottom friction vs Viscosity")

ax.set_title("a",loc="left",fontweight="bold")
ax1.set_title("b",loc="left",fontweight="bold")

plt.tight_layout()
plt.savefig(outpath+'plots/energy_budget_terms.png',dpi=200)
plt.close(fig)