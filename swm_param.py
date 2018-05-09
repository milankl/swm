## SET PARAMETERS
#
def set_param():
    """ seting up everything necessary for: parameters, grid, coriolis, forcing, initial conditions
        param = set_param() returns a parameter dictionary."""
    global param
    param = dict()

    ## parameters
    param['nx'] = 256               # number of grid points in x-direction
    param['ny'] = 256               # number of grid points in y-direction

    param['Lx'] = 3840e3            # basin width L [meters]
    param['Ly'] = 3840e3            # north-south basin extent [meters]

    param['g'] = 10.                # gravitational acceleration [ms**-2]
    param['H'] = 500.               # water depth [m]

    param['cfl'] = .9               # desired CFL-criterion, recommended 0.9, uses RK4
    param['Ndays'] = 5*365             # number of days to integrate

    param['dat_type'] = np.float32  # single/double precision use np.float32 or np.float64

    # initial conditions
    param['initial_conditions'] = 'rest'     # 'rest' or 'ncfile'
    param['init_run_id'] = 0                 # only for starting from ncfile
    param['init_interpolation'] = 0          # for 'ncfile': allow initial interpolation in case grids do not match

    # boundary conditions
    param['lbc'] = 0                         # no-slip: lbc=2, free-slip: lbc=0, 0<lbc<2 means partial-slip

    # bottom friction coefficient
    param['c_D'] = 1e-5

    # time stepping - do not change.
    param['scheme'] = 'RK4'

    # OUTPUT - of netcdf4, info_txt, parameters and scripts
    param['output'] = 1             # or 0 for no data storage
    param['output_dt'] = 6*3600    # every hours*3600 therefore in seconds
    param['outputpath'] = '/network/aopp/cirrus/pred/kloewer/swm_bf_cntrl/' # sets the path for netcdf output, else choose ''

    ## SET UP derived parameters
    set_grid()
    set_viscosity()
    set_coriolis()
    set_timestep()
    set_output()

    ## SET UP OPERATORS and FORCING
    set_grad_mat()      # set up the gradient matrices and make them globally available
    set_lapl_mat()      # set up harmonic and biharmonic diffusion (laplacians)
    set_interp_mat()    # set up the interpolation matrices and make them globally available
    set_arakawa_mat()   # set up the interpolation matrices for the Arakawa and Lamb scheme
    set_forcing()       # sets the wind forcing

    u,v,eta = initial_conditions()
    return u,v,eta

## grid parameters
def set_grid():
    """ The model is based on an Arakawa C-grid, with 4 staggered grids:

        T-grid: for eta, sits in the middle of a grid cell.
        u-grid: for u-velocities, sits in the middle of east&west edges
        v-grid: for v-velocities, sits in the middle of north&south edges
        q-grid: for vorticity, sits on corners of grid cells.

    """
    param['dx'] = param['Lx'] / param['nx'] # grid spacing in x-direction
    param['dy'] = param['Ly'] / param['ny'] # grid spacing in y-direction
    param['dA'] = param['dx']*param['dy']   # area of one grid cell
    param['lat_0'] = 30.                    # central latitude of beta-plane approximation [deg]
    param['max_dxdy'] = max(param['dx'],param['dy'])
    param['min_nxny'] = min(param['nx'],param['ny'])

    # feedback on grid
    param['a'] = 111194.   # 1 deg latitude in meters
    param['a_at_lat0'] = param['a']*np.cos(np.pi*param['lat_0']/180.)

    param['NT'] = param['nx']*param['ny']   # number of T-points (for h and tracers)
    param['Nu'] = (param['nx']-1)*param['ny']   # number of u-points
    param['Nv'] = param['nx']*(param['ny']-1)   # number of v-points
    param['Nq'] = (param['nx']+1)*(param['ny']+1) # number of q-points

    # grid vectors for T-points
    param['x_T'] = np.arange(param['dx']/2.,param['Lx'],param['dx'])
    param['y_T'] = np.arange(param['dy']/2.,param['Ly'],param['dy'])

    # grid vectors for u-points
    param['x_u'] = param['x_T'][:-1] + param['dx']/2.
    param['y_u'] = param['y_T']

    #grid vectors for v-points
    param['x_v'] = param['x_T']
    param['y_v'] = param['y_T'][:-1] + param['dy']/2.

    # grid vectors for q-points
    param['x_q'] = np.arange(0,param['Lx']+param['dx']/2.,param['dx'])
    param['y_q'] = np.arange(0,param['Ly']+param['dy']/2.,param['dy'])

## SET UP THE CORIOLIS MATRICES
def set_coriolis():
    """Sets up the coriolis parameter with beta-plane approximation as vector on the u-, v-, T- and q-grid."""

    global f_u,f_v,f_q,f_T

    # unpack dictionary for readability
    y_u, y_v = param['y_u'], param['y_v']
    y_T, y_q = param['y_T'], param['y_q']
    nx = param['nx']
    Ly = param['Ly']

    omega = 2*np.pi/(24.*3600.)     # Earth's angular frequency [s**-1]
    R = 6.371e6                     # Earth's radius [m]

    f_0 = 2*omega*np.sin(param['lat_0']*np.pi/180.)
    beta = 2*omega/R*np.cos(param['lat_0']*np.pi/180.)

    # store in dictionary
    param['f_0'] = f_0
    param['beta'] = beta

    # subtract the regions mid-y so that phi_0 corresponds to a central latitude
    yy_u = np.array([y_u - Ly/2.]*(nx-1)).T
    yy_v = np.array([y_v - Ly/2.]*nx).T
    yy_q = np.array([y_q - Ly/2.]*(nx+1)).T
    yy_T = np.array([y_T - Ly/2.]*nx).T

    # globally available coriolis parameters (only f_q is actually needed though)
    f_u = (f_0 + beta*yy_u.flatten()).astype(param['dat_type'])
    f_v = (f_0 + beta*yy_v.flatten()).astype(param['dat_type'])
    f_q = (f_0 + beta*yy_q.flatten()).astype(param['dat_type'])
    f_T = (f_0 + beta*yy_T.flatten()).astype(param['dat_type'])

## FORCING FIELDS
def set_forcing():
    """ Sets up the forcing field Fx of shape Nu and makes them globally available.
    This forcing is constant in time and x-direction, zero in y-direction, and varies only with y.
    Resembles trade winds and westerlies for a double gyre set up taken from Cooper and Zanna, 2015, Ocean Modelling.

    Note: the vector Fx excludes the 1/h-factor which is included in each model time step.
    """
    global Fx

    Lx,Ly = param['Lx'],param['Ly']     # for convenience
    param['rho'] = 1e3                  # density

    xx_u,yy_u = np.meshgrid(param['x_u'],param['y_u'])

    param['Fx0'] = 0.12      # was 0.12
    Fx = param['Fx0']*(np.cos(2*np.pi*(yy_u-Ly/2)/Ly) + 2*np.sin(np.pi*(yy_u - Ly/2)/Ly)) / param['rho']

    # from matrix to vector and set data type to have single/double precision
    Fx = Fx.flatten().astype(param['dat_type'])

## SET TIME STEPPING
def set_timestep():
    # shallow water gravity wave phase speed
    param['c_phase'] = np.sqrt(param['g']*param['H'])

    # the model timestep dt based on cfl stability criterion to resolve gravity waves
    # converting to integer, i.e. rounding down (floor)
    param['dt'] = (param['cfl']*min([param['dx'],param['dy']]) / param['c_phase']).astype(np.int64)

    # number of time steps to integrate
    param['Nt'] = np.ceil((param['Ndays'] * 3600. * 24.) / param['dt']).astype(np.int64)

## SET OUTPUT
def set_output():
    """ Creates folder for ouput. Initializes the nc-files etc. """

    param['path'] = path

    # output every n time steps
    # due to int the output time step might be a bit less than desired
    param['output_n'] = (param['output_dt']/param['dt']).astype(np.int64)
    param['output_tlen'] = np.ceil(param['Nt']/float(param['output_n'])).astype(np.int64)

    # hard coded fix for the case when Nt is a multiplicative of output_n
    if param['output_tlen']*param['output_n'] == param['Nt']:
        param['output_tlen'] += 1

    # set up next run folder run????
    if param['output']:
        cwd = os.getcwd()   # get current directory for switching back later
        if param['outputpath'] == '': # no specific outputpath provided
            try:
                os.chdir(param['path']+'/data')
            except:
                os.mkdir(param['path']+'/data')
                os.chdir(param['path']+'/data')
        else:
            os.chdir(param['outputpath'])

        all_runs = glob.glob('run????')
        if not all_runs:    # empty list
            param['run_id'] = 0
        else:               # add a new run id by taking the largest existing number and +1
            param['run_id'] = max([int(run[3:]) for run in all_runs])+1

        param['runfolder'] = 'run%04i' % param['run_id']
        os.mkdir(param['runfolder'])
        if param['outputpath'] == '':
            param['output_runpath'] = param['path']+'data/'+param['runfolder']
        else:
            param['output_runpath'] = param['outputpath']+param['runfolder']

        os.chdir(cwd)       # switch back to old directory

        # Save grid information in txt file
        output_txt_ini()
        str_tmp1 = (param['Lx']/1e3,param['Ly']/1e3)
        str_tmp2 = (param['Lx']/param['a_at_lat0'],param['Ly']/param['a'],param['lat_0'])
        str_tmp3 = (param['dx']/1e3,param['dy']/1e3)
        str_tmp4 = (param['dx']/param['a_at_lat0'],param['dy']/param['a'])

        output_txt('Domain is %ikm x %ikm (approx. %.1fdeg x %.1fdeg) centred at %.1fdegN' % (str_tmp1+str_tmp2))
        output_txt('with %i x %i grid points' % (param['nx'],param['ny']))
        output_txt('at resolution %.2fkm x %.2fkm (approx. %.2fdeg x %.2fdeg)' % (str_tmp3+str_tmp4))
        output_txt('')

        # Save all scripts as zipped file
        output_scripts()

## INITIALIZE PROGNOSTIC VARIABLES u,v,eta
def initial_conditions():
    """ Preallocates and sets the initial conditions for u,v,eta. """

    if param['initial_conditions'] == 'rest':
        u_0 = np.zeros(param['Nu']).astype(param['dat_type'])
        v_0 = np.zeros(param['Nv']).astype(param['dat_type'])
        eta_0 = (np.zeros(param['NT'])).astype(param['dat_type'])
        param['t0'] = 0

    elif param['initial_conditions'] == 'ncfile':
        if param['outputpath'] == '':
            initpath = param['path']+'data/run%04i' % param['init_run_id']
        else:
            initpath = param['outputpath']+'run%04i' % param['init_run_id']

        init_ncu = Dataset(initpath+'/u.nc')
        init_ncv = Dataset(initpath+'/v.nc')
        init_nceta = Dataset(initpath+'/eta.nc')

        u_0 = init_ncu['u'][-1,:,:]
        v_0 = init_ncv['v'][-1,:,:]
        eta_0 = init_nceta['eta'][-1,:,:]
        param['t0'] = init_nceta['t'][-1]
        output_txt('Starting from last state of run %04i' % param['init_run_id'])

        u_0 = u_0.flatten().astype(param['dat_type'])
        v_0 = v_0.flatten().astype(param['dat_type'])
        eta_0 = eta_0.flatten().astype(param['dat_type'])

    return u_0,v_0,eta_0

def set_viscosity():
    """ linear scaling of constant viscosity coefficients based on

        nu_A = 540 m**2/s at dx = 30km."""

    param['nu_A'] = 128*540./param['min_nxny']              # harmonic mixing coefficient
    param['nu_B'] = param['nu_A']*param['max_dxdy']**2      # biharmonic mixing coefficient
