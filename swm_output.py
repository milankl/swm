## OUTPUT FUNCTIONS
# PART1: STORE DATA in netCDF4 file (output__nc_ini,output_nc,output_nc_fin)
# PART2: STORE INFO in txt file (output_txt_ini, ...
# PART3: STORE PARAMETERS IN .NPY FILE

## STORE DATA
def output_nc_ini():
    """ Initialise the netCDF4 file."""

    param['output_j'] = 0   # output index

    # store files, dimensions and variables in dictionnaries
    ncu = dict()
    ncv = dict()
    nceta = dict()

    # creating the netcdf files
    ncformat = 'NETCDF4'
    ncu['file'] = Dataset(param['output_runpath']+'/u.nc','w',format=ncformat)
    ncv['file'] = Dataset(param['output_runpath']+'/v.nc','w',format=ncformat)
    nceta['file'] = Dataset(param['output_runpath']+'/eta.nc','w',format=ncformat)

    # write general attributes
    for ncfile in [ncu,ncv,nceta]:
        ncfile['file'].history = 'Created ' + tictoc.ctime(tictoc.time())
        ncfile['file'].description = 'Data from: Shallow-water model in double gyre configuration.'
        ncfile['file'].details = 'Cartesian coordinates, beta-plane approximation, Arakawa C-grid'

        # all param ints floats and strings as global attribute
        for key in param.keys():
            if (type(param[key]) is int) or (type(param[key]) is float) or (type(param[key]) is str):
                ncfile['file'].setncattr(key,param[key])

    # create dimensions
    ncu['xdim'] = ncu['file'].createDimension('x',param['nx']-1)
    ncu['ydim'] = ncu['file'].createDimension('y',param['ny'])
    #ncu['tdim'] = ncu['file'].createDimension('t',param['output_tlen'])
    ncu['tdim'] = ncu['file'].createDimension('t',None)


    ncv['xdim'] = ncv['file'].createDimension('x',param['nx'])
    ncv['ydim'] = ncv['file'].createDimension('y',param['ny']-1)
    ncv['tdim'] = ncv['file'].createDimension('t',None)

    nceta['xdim'] = nceta['file'].createDimension('x',param['nx'])
    nceta['ydim'] = nceta['file'].createDimension('y',param['ny'])
    nceta['tdim'] = nceta['file'].createDimension('t',None)

    # create variables
    p = 'f4' # 32-bit precision storing, or f8 for 64bit
    for ncfile,var in zip([ncu,ncv,nceta],['u','v','eta']):
        # store time as integers as measured in seconds and gets large
        ncfile['t'] = ncfile['file'].createVariable('t','i8',('t',),zlib=True,fletcher32=True)
        ncfile['x'] = ncfile['file'].createVariable('x','f8',('x',),zlib=True,fletcher32=True)
        ncfile['y'] = ncfile['file'].createVariable('y','f8',('y',),zlib=True,fletcher32=True)
        ncfile[var] = ncfile['file'].createVariable(var,p,('t','y','x'),zlib=True,fletcher32=True)

    # write units
    for ncfile in [ncu,ncv,nceta]:
        ncfile['t'].units = 's'
        ncfile['t'].long_name = 'time'
        ncfile['x'].units = 'm'
        ncfile['x'].long_name = 'x'
        ncfile['y'].units = 'm'
        ncfile['y'].long_name = 'y'

    ncu['u'].units = 'm/s'
    ncv['v'].units = 'm/s'
    nceta['eta'].units = 'm'

    # write dimensions
    for ncfile,var in zip([ncu,ncv,nceta],['u','v','T']):
        ncfile['x'][:] = param['x_'+var]
        ncfile['y'][:] = param['y_'+var]

    # make globally available
    global ncfiles
    ncfiles = [ncu,ncv,nceta]

    output_txt('Output will be stored in '+param['outputpath']+param['runfolder']+' every %i hours.' % (param['output_dt']/3600.))


def output_nc(u,v,eta,t):
    """ Writes u,v,eta fields on every nth time step """
    # output index j
    j = param['output_j']   # for convenience

    for ncfile in ncfiles:
        ncfile['t'][j] = t

    #TODO issue, use unlimited time dimension or not?
    ncfiles[0]['u'][j,:,:] = u2mat(u)
    ncfiles[1]['v'][j,:,:] = v2mat(v)
    ncfiles[2]['eta'][j,:,:] = h2mat(eta)

    param['output_j'] += 1

def output_nc_fin():
    """ Finalise the output netCDF4 file."""

    for ncfile in ncfiles:
        ncfile['file'].close()

    output_txt('All output written in '+param['runfolder']+'.')

## STORE INFO in TXT FILE
def readable_secs(secs):
    """ Returns a human readable string representing seconds in terms of days, hours, minutes, seconds. """

    days = np.floor(secs/3600/24)
    hours = np.floor((secs/3600) % 24)
    minutes = np.floor((secs/60) % 60)
    seconds = np.floor(secs%3600%60)

    if days > 0:
        return ("%id, %ih" % (days,hours))
    elif hours > 0:
        return ("%ih, %imin" % (hours,minutes))
    elif minutes > 0:
        return ("%imin, %is" % (minutes,seconds))
    else:
        return ("%.2fs" % secs)

def duration_est(tic):
    """ Saves an estimate for the total time the model integration will take in the output txt file. """
    time_togo = (tictoc.time()-tic) / (i+1) * param['Nt']
    str1 = 'Model integration will take approximately '+readable_secs(time_togo)+', '
    print(str1)

    if param['output']:
        str2 = 'and is hopefully done on '+tictoc.asctime(tictoc.localtime(tic + time_togo))
        output_txt(str1+str2)
        print(str2)

def output_txt_ini():
    """ Initialise the output txt file for information about the run."""
    if param['output']:
        param['output_txtfile'] = open(param['output_runpath']+'/info.txt','w')
        s = ('Shallow water model run %i initialised on ' % param['run_id'])+tictoc.asctime()+'\n'
        param['output_txtfile'].write(s)

def output_scripts():
    """Save all model scripts into a zip file."""
    if param['output']:
        zf = zipfile.ZipFile(param['output_runpath']+'/scripts.zip','w')
        all_scripts = glob.glob('swm_*.py')
        [zf.write(script) for script in all_scripts]
        zf.close()
        output_txt('All model scripts stored in a zipped file.')

def output_txt(s,end='\n'):
    """ Write into the output txt file."""
    if param['output']:
        param['output_txtfile'].write(s+end)
        param['output_txtfile'].flush()

def output_txt_fin():
    """ Finalise the output txt file."""
    if param['output']:
        param['output_txtfile'].close()

## STORE PARAMETERS
def output_param():
    """ Stores the param dictionary in a .npy file """
    if param['output']:
        # filter out 'output_txtfile' as this is a unsaveable textwrapper
        dict_tmp = {key:param[key] for key in param.keys() if key != 'output_txtfile'}
        np.save(param['output_runpath']+'/param.npy',dict_tmp)

        # store also as a more readable .txt file for quick access on the parameters
        param_txtfile = open(param['output_runpath']+'/param.txt','w')
        for key in dict_tmp.keys():
            if not(key in ['x_T','y_T','x_u','y_u','x_v','y_v','x_q','y_q']):
                param_txtfile.write(key + 2*'\t' + str(dict_tmp[key]) + '\n')

        param_txtfile.close()
        output_txt('Param dictionary stored as txt and zip.\n')
