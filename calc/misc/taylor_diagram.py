## TAYLOR DIAGRAM PLOTTING FUNCTION:

def taylor_diag(ref,mod,norm=False,units='',label='mod',lsize=11,marker='^',color='r',msize=5,alpha=1,title='Taylor Diagram',lim=None,ax=None,add=False,return_stats=False):
    """
    plot the taylor diagram showing the skill of the modeled field 'mod' to simulate the observed field 'ref', note that this   excludes all information about biases in the simulated field!
    
    INPUT:
        ref: array, reference/observed field
        mod: array, simulated field, must have same dimensions as ref
        norm: True/False, if norm=True normalize the standard deviations and RMS so that fields of different variables (different units) can be plotted into the same diagram
        units: string for axes' labels, only necessary if norm=False
        weights: array of same dimension as ref, if specified, weights will be used to calculate mean, std and RMS
        lim: upper axes limit, should be greater than std of ref! if not specified it will be computed automatically
            CAREFUL: right now, the limit will be redefined within the function so the actual limit will never be lim...
        ax: axis that should be plotted onto, if not specified a new plot will be set up
        add: True/False, if True, a single point will be added to the existing Taylor Diagram that must be specified with the axis handle ax!
    
    OUTPUT:
        returns the Taylor diagram axis handle (always) and figure handle (if ax=None)
    """
    # transparency of the axes:
    alp = .5
    
    if (ax is None and add):
        print('pass handle for an existing axis in ax to add the point!')
        return None
        
    ref_std = ref.std()
    mod_std = mod.std()
    R = np.corrcoef(ref,mod)[0,1]
    MSD = ((ref-mod)**2).mean()
    
    if norm:
        # normalize by std of the reference field
        mod_std /= ref_std; MSD /= ref_std**2
        ref_std /= ref_std
        units = ''
        
    if add==False:
        # SET UP PLOT (if necessary):
        if ax is None:
            fig = plt.figure(figsize=(6,6))
            # set up the Taylor diagram:
            ax = plt.subplot(111)
    
        
        # automatically choose a LIMIT for the plot and define the std circle lengths:
        decpl = np.ceil(np.log10(ref_std))
        if lim is None:
            lim = np.round(4/3*ref_std,int(-decpl)+1)
            lim = 1.1
        # ncirc = lim/(10**(decpl-1))
        # clens = np.linspace(0,lim,num=ncirc+1)
        
        # set limits, aspect ratio, labels and axes
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        ax.set_ylim([0,1.05*lim])
        ax.set_xlim([0,1.05*lim])
        ax.set_aspect(1)
        ax.set_ylabel('Standard Deviation ' + units)
        ax.set_xlabel('Standard Deviation ' + units)
        clens = [1.2,0.8,0.4]
        lim = clens[0]
        ax.set_ylim([0,1.05*lim])
        ax.set_xlim([0,1.05*lim])

        # STANDARD DEVIATION
        std_circ_angles = np.arange(np.pi/2.,-.001,-.001)
        for (CL,cc) in zip(clens,range(len(clens))):
            x_circ_o, y_circ_o = polar2z(CL,std_circ_angles)
            if cc==0:
                ax.plot(x_circ_o,y_circ_o,color='g',alpha=alp)
            else:    
                ax.plot(x_circ_o,y_circ_o,color='k',alpha=alp)
            
        x_circ,y_circ = polar2z(ref_std,std_circ_angles)
        ax.plot(x_circ,y_circ,color='k',linestyle='dashed',alpha=alp)
        
        # CORRELATION
        corrvals = np.arange(0.0,1.01,.1)
        corrvals_tick = np.concatenate((np.arange(0.0,.91,.05),np.arange(.91,1.,.01)))
        corr_ex_ticks = np.array([0.95,0.99])

        x_corr, y_corr = polar2z(lim,np.arccos(corrvals))
        x_corr_tex,y_corr_tex = polar2z(1.04*lim,np.arccos(corrvals))
        x_corr_tex_ex, y_corr_tex_ex = polar2z(1.04*lim,np.arccos(corr_ex_ticks))
        x_tick1, y_tick1 = polar2z(.99*lim,np.arccos(corrvals_tick))
        x_tick2, y_tick2 = polar2z(lim,np.arccos(corrvals_tick))
        
        [ax.plot([0,x_corr[ii]],[0,y_corr[ii]],color='g',alpha=alp) for ii in range(len(corrvals))]
        [ax.plot([x_tick1[jj],x_tick2[jj]],[y_tick1[jj],y_tick2[jj]],color='g',alpha=alp) for jj in range(len(corrvals_tick))]
        [ax.text(x_corr_tex[kk],y_corr_tex[kk],str(CV),color='g',alpha=alp,rotation=np.rad2deg(np.arccos(CV)),va='center',ha='center') for (CV,kk) in zip(corrvals,range(len(corrvals)))]
        [ax.text(x_corr_tex_ex[ll],y_corr_tex_ex[ll],str(CEXT),color='g',alpha=alp,rotation=np.rad2deg(np.arccos(CEXT)),va='center',ha='center') for (CEXT,ll) in zip(corr_ex_ticks,range(len(corr_ex_ticks)))]
        ax.text(.78*lim,.78*lim,'Correlation',rotation=315,fontsize=15,ha='center',va='center',color='g',alpha=alp)
        
        # ROOT MEAN SQUARE DEVIATION
        for CLE in clens:
            if CLE!=0:
                lamb = np.arange(0,np.pi,.001)
                x_RMS,y_RMS = polar2z(CLE,lamb)
                x_RMS += ref_std
                inlim = np.where(np.logical_and(z2polar(x_RMS,y_RMS)[0]<lim,x_RMS>0))
                lamb = lamb[inlim]
                y_RMS = y_RMS[inlim]
                x_RMS = x_RMS[inlim]
                RMSloc = int(len(x_RMS)/2)
            
            ax.plot(x_RMS,y_RMS,color='y',alpha=alp,linestyle='solid')
            ax.text(x_RMS[RMSloc],y_RMS[RMSloc],'%.1f' % CLE,color='y',alpha=alp,rotation=np.rad2deg(lamb[RMSloc])+270,ha='center',va='center',bbox=dict(facecolor='white', edgecolor='none',pad=0.1,boxstyle='round'))
        
        # show the reference std by a circle marked 'obs'
        ax.plot(ref_std,0,'o',color='C2',clip_on=False)
        #ax.text(ref_std,-.03,'ref',ha='center',va='top')
        AXTE = [AXT for AXT in ax.get_xticks() if not np.allclose(AXT,ref_std,rtol=.05)]
        ax.set_xticks(AXTE)
        ax.set_yticks(ax.get_yticks())
        
    # plot the point representing the model fields statistics:
    x_mod,y_mod = polar2z(mod_std,np.arccos(R))
    ax.plot(x_mod,y_mod,marker,color=color,markeredgecolor=color,markersize=msize,alpha=alpha,clip_on=False,zorder=5)
    ax.text(x_mod,y_mod,' ' + label,fontsize=lsize,color=color,zorder=10)
    ax.set_title(title,fontsize=16)
    
    if return_stats:
        try:
            return fig,ax,(ref_std,mod_std,R,MSD)
        except NameError:
            return ax,(ref_std,mod_std,R,MSD)
    else:
        try:
            return fig,ax
        except NameError:
            return ax
            
    

def polar2z(r,theta):
    eutf = r * np.exp( 1j * theta )
    return np.real(eutf),np.imag(eutf)

def z2polar(x,y):
    z = x + 1j*y
    return np.abs(z),np.angle(z)

def grid_interpolation(u,v,h,param_old,param_new):
    from scipy.interpolate import RegularGridInterpolator as RIG
    
    # from vector to matrix
    u = u2mat(u,param_old)
    v = v2mat(v,param_old)
    h = h2mat(h,param_old)
    
    # padding following the boundary conditions (kinematic + no-slip)
    u = np.pad(u,((1,1),(1,1)),'constant',constant_values=((0,0),(0,0)))    # kinematic bc
    v = np.pad(v,((1,1),(1,1)),'constant',constant_values=((0,0),(0,0)))
    h = np.pad(h,1,'edge')                  # no gradients of h across boundaries
    
    # padding also for x,y grid
    x_u_old = np.hstack((0,param_old['x_u'],param_old['Lx']))
    y_u_old = np.hstack((0,param_old['y_u'],param_old['Ly']))
    
    x_v_old = np.hstack((0,param_old['x_v'],param_old['Lx']))
    y_v_old = np.hstack((0,param_old['y_v'],param_old['Ly']))
    
    x_T_old = np.hstack((0,param_old['x_T'],param_old['Lx']))
    y_T_old = np.hstack((0,param_old['y_T'],param_old['Ly']))
    
    # setting up the interpolation functions
    RIG_u = RIG((y_u_old,x_u_old),u)
    RIG_v = RIG((y_v_old,x_v_old),v)
    RIG_h = RIG((y_T_old,x_T_old),h)
    
    xx_T,yy_T = np.meshgrid(param_new['x_T'],param_new['y_T'])
    xx_u,yy_u = np.meshgrid(param_new['x_u'],param_new['y_u'])
    xx_v,yy_v = np.meshgrid(param_new['x_v'],param_new['y_v'])
    
    # receive the interpolated initial conditions
    u_new = RIG_u((yy_u.flatten(),xx_u.flatten()))
    v_new = RIG_v((yy_v.flatten(),xx_v.flatten()))
    h_new = RIG_h((yy_T.flatten(),xx_T.flatten()))
    
    return u_new,v_new,h_new
 
def grid_interpolation_T(T,param_old,param_new):
    from scipy.interpolate import RegularGridInterpolator as RIG
    
    # from vector to matrix
    T = h2mat(T,param_old)
    T = np.pad(T,1,'edge')                  # no gradients of T across boundaries
    
    # padding also for x,y grid
    x_T_old = np.hstack((0,param_old['x_T'],param_old['Lx']))
    y_T_old = np.hstack((0,param_old['y_T'],param_old['Ly']))
    
    # setting up the interpolation functions
    RIG_T = RIG((y_T_old,x_T_old),T)    
    xx_T,yy_T = np.meshgrid(param_new['x_T'],param_new['y_T'])
    
    # receive the interpolated values
    T_new = RIG_T((yy_T.flatten(),xx_T.flatten()))
    
    return T_new
    
def h2mat(h,param):
    return h.reshape((param['ny'],param['nx']))

def u2mat(u,param):
    return u.reshape((param['ny'],param['nx']-1))

def v2mat(v,param):
    return v.reshape((param['ny']-1,param['nx']))