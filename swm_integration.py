## TIME INTEGRATION
def time_integration(u,v,h):
    
    tic = tictoc.time()     # measure time
    global dt
    dt = param['dt']        # for convenience
    t = param['t0']         # initial time
    feedback_ini(u,v,h,t)   # output

    ## RUNGE KUTTA 3rd ORDER or 4th ORDER
    if param['scheme'][:2] == 'RK':
        if param['scheme'] == 'RK4':    # rename either 3rd order or 4th scheme RKn
            RKn = RK4
        elif param['scheme'] == 'RK3':
            RKn = RK3
    
        # same loop for RK4 and for RK3 due to renaming above
        global i
        for i in range(param['Nt']):
            duvh = RKn(u,v,h)
            u = u + dt*duvh[0]
            v = v + dt*duvh[1]
            h = h + dt*duvh[2]
            t += dt
            
            feedback(u,v,h,t,tic)

    ## ADAMS-BASHFORTH
    else:
        NAB = int(param['scheme'][2:])      # order of adams-bashforth multistep scheme
        swap = [NAB-1]+list(range(0,NAB-1)) # indices to swap the last rhs after each time step
        ABb = ABcoefficients(NAB)           # ADAMS-BASHFORTH coefficients
        
        # preallocate the multistep matrices (right-hand sides for the last NABth time steps)
        urhs = np.empty((param['Nu'],NAB)).astype(param['dat_type'])
        vrhs = np.empty((param['Nv'],NAB)).astype(param['dat_type'])
        hrhs = np.empty((param['NT'],NAB)).astype(param['dat_type'])
        
        for i in range(param['Nt']):
            # increases the AB order until the desired order is reached
            abcolumn = min(i,NAB-1) 

            urhs[:,0],vrhs[:,0],hrhs[:,0] = rhs(u,v,h)  # evaluation of the rhs
            u = u + dt*urhs.dot(ABb[:,abcolumn])        # update u,v,h
            v = v + dt*vrhs.dot(ABb[:,abcolumn])
            h = h + dt*hrhs.dot(ABb[:,abcolumn])
            t += dt

            urhs = urhs[:,swap]     # swap multistep matrices
            vrhs = vrhs[:,swap]     # e.g. 0->1, 1->2, 2->3, 3->0
            hrhs = hrhs[:,swap]
            
            feedback(u,v,h,t,tic)
    
    print(('Integration done in '+readable_secs(tictoc.time() - tic)+' on '+tictoc.asctime()))
    output_txt(('\nTime integration done in '+readable_secs(tictoc.time() - tic)+' on '+tictoc.asctime()))

    # finalising output
    if param['output']:
        output_nc_fin()         # finalise nc file
        output_txt_fin()        # finalise info txt file
        
    return u,v,h

### TIME STEPPING SCHEMES
def RK4(u,v,h):
    """ Computes the right-hand side using RUNGE KUTTA 4th order scheme. 
    u,v,h are coupled in every of the 4 sub-time steps of the RK4 scheme."""
    k1 = rhs(u,v,h)
    k2 = rhs(u + dt/2.*k1[0],v + dt/2.*k1[1],h + dt/2.*k1[2])
    k3 = rhs(u + dt/2.*k2[0],v + dt/2.*k2[1],h + dt/2.*k2[2])
    k4 = rhs(u + dt*k3[0],v + dt*k3[1],h + dt*k3[2])
    
    du = (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]) / 6.
    dv = (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]) / 6.
    dh = (k1[2] + 2*k2[2] + 2*k3[2] + k4[2]) / 6.
    
    return du,dv,dh
    
def RK3(u,v,h):
    """ Computes the right-hand side using RUNGE KUTTA 3rd order scheme. 
    u,v,h are coupled in every of the 3 sub-time steps of the RK3 scheme."""
    k1 = rhs(u,v,h)
    k2 = rhs(u + dt/2.*k1[0],v + dt/2.*k1[1],h + dt/2.*k1[2])
    k3 = rhs(u + dt*(-k1[0]+2*k2[0]),v + dt*(-k1[1]+2*k2[1]),h + dt*(-k1[2]+2*k2[2]))
    
    du = (k1[0] + 4*k2[0] + k3[0]) / 6.
    dv = (k1[1] + 4*k2[1] + k3[1]) / 6.
    dh = (k1[2] + 4*k2[2] + k3[2]) / 6.
    
    return du,dv,dh
    
def ABcoefficients(N):
    """ Returns the Adams-Bashforth coefficients up to order N <= 5. """
    ABb = np.array([[ 1.        ,  1.5       ,     23./12.,     55./24.,  1901./720.],\
                    [ 0.        , -0.5       ,      -4./3.,    -59./24., -1387./360.],\
                    [ 0.        ,  0.        ,      5./12.,     37./24.,    109./30.],\
                    [ 0.        ,  0.        ,  0.        ,       -3/8.,  -637./360.],\
                    [ 0.        ,  0.        ,  0.        ,  0.        ,   251./720.]])
                        
    return ABb[:N,:N].astype(param['dat_type'])

## FEEDBACK ON INTEGRATION
def feedback_ini(u,v,h,t):
    if param['output']:
        output_nc_ini()
        output_nc(u,v,h,t)  # store initial conditions
        output_param()      # store the param dictionnary
        
        # Store information in txt file
        output_txt('Integrating %.1f days with dt=%.2f min in %i time steps' % (param['Ndays'],dt/60.,param['Nt']))    
        output_txt('and eddy viscosities A = %4.1f, B = %4.1f' % (param['A'],param['B']))
        output_txt('Time integration scheme is '+param['scheme']+' with CFL = %.2f on %i threads.' % (param['cfl'],param['num_threads']))
        output_txt('')
        output_txt('Starting shallow water model on '+tictoc.asctime())
        print(('Starting shallow water model run %i on ' % param['run_id'])+tictoc.asctime())
    else:
        print('Starting shallow water model on '+tictoc.asctime())

def feedback(u,v,h,t,tic):
    if (i+1) % param['output_n'] == 0:
        if param['output']:     # storing u,v,h as netCDF4
            output_nc(u,v,h,t)
    
    # feedback on progress every integer % step.
    if ((i+1)/param['Nt']*100 % 1) < (i/param['Nt']*100 % 1):
        progress = str(int((i+1)/param['Nt']*100.))+'%'
        print(progress, end='\r')
        if i > 100:
            output_txt(progress,'\n')

    if i == 100: 
        # estimate total time for integration after 100 time steps.
        duration_est(tic)
