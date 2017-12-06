## TIME INTEGRATION
def time_integration(u,v,eta,e):

    tic = tictoc.time()     # measure time
    global dt
    dt = param['dt']        # for convenience
    t = param['t0']         # initial time
    feedback_ini(u,v,eta,e,t)   # output

    global i
    for i in range(param['Nt']):
        duvetae = RK4(u,v,eta,e)    # only RK4 allowed - other explicit schemes are inferior
        u += dt*duvetae[0]
        v += dt*duvetae[1]
        eta += dt*duvetae[2]
        e += dt*duvetae[3]
        t += dt

        feedback(u,v,eta,e,t,tic)

    print(('Integration done in '+readable_secs(tictoc.time() - tic)+' on '+tictoc.asctime()))
    output_txt(('\nTime integration done in '+readable_secs(tictoc.time() - tic)+' on '+tictoc.asctime()))

    # finalising output
    if param['output']:
        output_nc_fin()         # finalise nc file
        output_txt_fin()        # finalise info txt file

    return u,v,eta,e

### TIME STEPPING SCHEMES
def RK4(u,v,eta,e):
    """ Computes the right-hand side using RUNGE KUTTA 4th order scheme.
    u,v,h are coupled in every of the 4 sub-time steps of the RK4 scheme."""
    k1 = rhs(u,v,eta,e)
    k2 = rhs(u + dt/2.*k1[0],v + dt/2.*k1[1],eta + dt/2.*k1[2],e + dt/2.*k1[3])
    k3 = rhs(u + dt/2.*k2[0],v + dt/2.*k2[1],eta + dt/2.*k2[2],e + dt/2.*k2[3])
    k4 = rhs(u + dt*k3[0],v + dt*k3[1],eta + dt*k3[2],e + dt*k3[3])

    du = (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]) / 6.
    dv = (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]) / 6.
    deta = (k1[2] + 2*k2[2] + 2*k3[2] + k4[2]) / 6.
    de = (k1[3] + 2*k2[3] + 2*k3[3] + k4[3]) / 6.

    return du,dv,deta,de

## FEEDBACK ON INTEGRATION
def feedback_ini(u,v,eta,e,t):
    if param['output']:
        output_nc_ini()
        output_nc(u,v,eta,e,t)  # store initial conditions
        output_param()      # store the param dictionnary

        # Store information in txt file
        output_txt('Integrating %.1f days with dt=%.2f min in %i time steps' % (param['Ndays'],dt/60.,param['Nt']))
        output_txt('Time integration scheme is RK4 with CFL = %.2f' % param['cfl'])
        output_txt('')
        output_txt('Starting shallow water model on '+tictoc.asctime())
        print(('Starting shallow water model run %i on ' % param['run_id'])+tictoc.asctime())
    else:
        print('Starting shallow water model on '+tictoc.asctime())

def feedback(u,v,eta,e,t,tic):
    if (i+1) % param['output_n'] == 0:
        if param['output']:     # storing u,v,h as netCDF4
            output_nc(u,v,eta,e,t)

    # feedback on progress every 5% step.
    if ((i+1)/param['Nt']*100 % 5) < (i/param['Nt']*100 % 5):
        progress = str(int((i+1)/param['Nt']*100.))+'%'
        print(progress, end='\r')
        if i > 100:
            output_txt(progress,'\n')

    if i == 100:
        # estimate total time for integration after 100 time steps.
        duration_est(tic)
