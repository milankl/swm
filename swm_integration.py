## TIME INTEGRATION
def time_integration(u,v,eta):

    tic = tictoc.time()     # measure time
    global dt
    dt = param['dt']        # for convenience
    t = param['t0']         # initial time
    feedback_ini(u,v,eta,t)   # output

    ## RUNGE KUTTA 4th ORDER
    rk_a = np.array([1/6.,1/3.,1/3.,1/6.])
    rk_b = np.array([0.5,0.5,1.])

    # can't trigger deep copy through [:] use .copy() instead
    u0,v0,eta0 = u.copy(),v.copy(),eta.copy()
    u1,v1,eta1 = u.copy(),v.copy(),eta.copy()

    global i    # iteration index
    for i in range(param['Nt']):

        # trigger deep copy through [:]
        u1[:],v1[:],eta1[:] = u,v,eta

        for rki in range(4):
            du,dv,deta = rhs(u1,v1,eta1)

            if rki < 3: # RHS update for the next RK-step
                u1 = u + rk_b[rki]*dt*du
                v1 = v + rk_b[rki]*dt*dv
                eta1 = eta + rk_b[rki]*dt*deta

            # Summing all the RHS on the go
            u0 += rk_a[rki]*dt*du
            v0 += rk_a[rki]*dt*dv
            eta0 += rk_a[rki]*dt*deta

        # deep copy through [:]
        u[:],v[:],eta[:] = u0,v0,eta0

        t += dt
        feedback(u,v,eta,t,tic)

    print(('Integration done in '+readable_secs(tictoc.time() - tic)+' on '+tictoc.asctime()))
    output_txt(('\nTime integration done in '+readable_secs(tictoc.time() - tic)+' on '+tictoc.asctime()))

    # finalising output
    if param['output']:
        output_nc_fin()         # finalise nc file
        output_txt_fin()        # finalise info txt file

    return u,v,eta

## FEEDBACK ON INTEGRATION
def feedback_ini(u,v,eta,t):
    if param['output']:
        output_nc_ini()
        output_nc(u,v,eta,t)  # store initial conditions
        output_param()      # store the param dictionnary

        # Store information in txt file
        output_txt('Integrating %.1f days with dt=%.2f min in %i time steps' % (param['Ndays'],dt/60.,param['Nt']))
        output_txt('Time integration scheme is '+param['scheme']+' with CFL = %.2f' % param['cfl'])
        output_txt('')
        output_txt('Starting shallow water model on '+tictoc.asctime())
        print(('Starting shallow water model run %i on ' % param['run_id'])+tictoc.asctime())
    else:
        print('Starting shallow water model on '+tictoc.asctime())

def feedback(u,v,eta,t,tic):
    if (i+1) % param['output_n'] == 0:
        if param['output']:     # storing u,v,h as netCDF4
            output_nc(u,v,eta,t)

    # feedback on progress every 5% step.
    if ((i+1)/param['Nt']*100 % 5) < (i/param['Nt']*100 % 5):
        progress = str(int((i+1)/param['Nt']*100.))+'%'
        print(progress, end='\r')
        if i > 100:
            output_txt(progress,'\n')

    if i == 100:
        # estimate total time for integration after 100 time steps.
        duration_est(tic)
