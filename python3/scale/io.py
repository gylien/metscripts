import numpy as np
import numpy.ma as ma
import ncphysio
import os.path
import datetime as dt
from netCDF4 import Dataset


__all__ = ['scale_dimlist', 'scale_dimlist_g', 'scale_file_suffix', 'scale_time_0',
           'scale_open', 'scale_close', 'scale_gettime', 'scale_puttime',
           'scale_read', 'scale_write', 'ScaleIO']


scale_dimlist = [
['time', 'time1'],
['nv'],
['z', 'zh', 'lz', 'lzh', 'uz', 'uzh', 'CZ', 'FZ', 'FDZ', 'LCZ', 'LFZ', 'UCZ', 'UFZ'],
['y', 'yh', 'CY', 'FY', 'FDY', 'CYG', 'FYG'],
['x', 'xh', 'CX', 'FX', 'FDX', 'CXG', 'FXG']
]
scale_dimlist_g = [
[],
[],
[],
['y', 'yh'],
['x', 'xh']
]
scale_file_suffix = '.pe{:06d}.nc'
scale_time_0 = dt.datetime(2013, 1, 1, 0, 0, 0)


def scale_open(basename, mode='r', scale_dimdef=None):
    """
    Open a set of split SCALE files and return the definition of 
    global and subdomain dimensions.

    Parameters
    ----------
    basename : string
        Split SCALE file basename. Path can be included.
    mode : string, optional
        File I/O mode: `r` for read and `r+` for read/write. Default: `r`
    scale_dimdef : array of array, optional
        List of dimensions in the SCALE files. Default: `scale_dimlist`

    Returns
    -------
    nproc : integer
        Number of split SCALE files
    rootgrps : array of netcdf4-python Dataset instance
        Array of netcdf4-python Dataset instances of the split SCALE files.
    scale_dimdef : dictionary
        Summary of dimensions in the split SCALE files
        scale_dimdef['len'] : dictionary
            Lengths of local dimensions in the split files
        scale_dimdef['len_g'] : dictionary
            Lengths of global dimensions
        scale_dimdef['coor_g'] : dictionary
            Coordinates of global dimensions
        scale_dimdef['start'] : dictionary
            Start indices of global dimensions
    """
    ncfile = basename + scale_file_suffix.format(0)
    if not os.path.isfile(ncfile):
        if mode == 'r':
            raise IOError("File does not exist... basename = '" + basename + "'")
        elif scale_dimdef is None:
            raise IOError("File does not exist... basename = '" + basename + "'")
        else:
#            scale_create_new(basename, scale_dimdef)
            raise IOError('Scale_create has not been supported yet...')

    rootgrps = []
    sub_ip = {}
    sub_idx = {}
    sub_var = {}
    ip = 0
    while True:
        ncfile = basename + scale_file_suffix.format(ip)
        if not os.path.isfile(ncfile):
            break

        rootgrps.append(Dataset(ncfile, mode))
        for idiml in scale_dimlist_g:
            for idim in idiml:
                if ip == 0:
                    sub_ip[idim] = [ip] * len(rootgrps[ip].dimensions[idim])
                    sub_idx[idim] = list(range(len(rootgrps[ip].dimensions[idim])))
                    sub_var[idim] = rootgrps[ip].variables[idim][:]
                else:
                    sub_ip[idim] += [ip] * len(rootgrps[ip].dimensions[idim])
                    sub_idx[idim] += list(range(len(rootgrps[ip].dimensions[idim])))
                    sub_var[idim] = np.append(sub_var[idim], rootgrps[ip].variables[idim][:])
        ip += 1
    nproc = ip

#    dimlen = ncphysio.ncphys_read_dimlen(rootgrps[0], dimlist=scale_dimlist)
    dimlen = {}
    for idiml in scale_dimlist:
        for idim in idiml:
            dimlen[idim] = [None] * nproc
            for ip in range(nproc): 
                if idim in rootgrps[ip].dimensions:
                    dimlen[idim][ip] = len(rootgrps[ip].dimensions[idim])
                else:
                    dimlen[idim][ip] = None

    dimcoor_g = {}
    dimlen_g = {}
    dimstart = {}
    for idiml in scale_dimlist_g:
        for idim in idiml:
            dimcoor_g[idim], indices = np.unique(sub_var[idim], return_inverse=True)
            dimlen_g[idim] = len(dimcoor_g[idim])
            dimstart[idim] = [None] * nproc
            for i, ip in enumerate(sub_ip[idim]):
                if dimstart[idim][ip] is None:
                    dimstart[idim][ip] = indices[i] - sub_idx[idim][i]
                elif dimstart[idim][ip] != indices[i] - sub_idx[idim][i]:
                    raise ValueError('Subdomains are not consistent.')

    scale_dimdef = {'len': dimlen, 'len_g': dimlen_g, 'coor_g': dimcoor_g, 'start': dimstart}
    return nproc, rootgrps, scale_dimdef


def scale_close(rootgrps):
    """
    Close a set of split SCALE files.

    Parameters
    ----------
    rootgrps : array of netcdf4-python Dataset instance
        Array of netcdf4-python Dataset instances of the split SCALE files.
    """
    for irg in rootgrps:
        irg.close()


def scale_gettime(scale_time):
    """
    Convert SCALE model time to python datetime.

    Parameters
    ----------
    scale_time : float
        Time in SCALE files

    Returns
    -------
    time : <datetime.datetime> class
        Time in <datetime.datetime> class
    """
    return scale_time_0 + dt.timedelta(seconds=scale_time)


def scale_puttime(time):
    """
    Convert python datetime to scale model time.

    Parameters
    ----------
    time : <datetime.datetime> class
        Time in <datetime.datetime> class

    Returns
    -------
    scale_time : float
        Time in SCALE files
    """
    return (time - scale_time_0).total_seconds()


def scale_read(nproc, rootgrps, scale_dimdef, varname, t=None):
    """
    Read a variable from a set of split SCALE files.

    Parameters
    ----------
    nproc : integer
        Number of split SCALE files
    rootgrps : array of netcdf4-python Dataset instance
        Array of netcdf4-python Dataset instances of the split SCALE files.
    scale_dimdef : dictionary
        Summary of dimensions in the split SCALE files
    varname : string
        The variable name.
    t : int or <datetime.datetime> class or None, optional
        Time to read. None for all times. Defalut: None

    Returns
    -------
    vardim : dictionary
        Dimensions of the return variable data.
    vardata : ndarray or masked_array
        Variable data in a ndarray or masked_array (if the variable has the 
        `_FillValue` attribute).
    """
    if t is None:
        it = 'all'
        time = 'all'
    elif type(t) is int:
        it = t
        time = None
    elif type(t) is dt.datetime:
        it = None
        time = scale_puttime(t)
    else:
        raise ValueError("The type of 't' should be either 'int' or 'datetime.datetime' or 'None'.")

    vardim, vardata_0 = ncphysio.ncphys_read(rootgrps[0], varname, dimlist=scale_dimlist, time=time, it=it)
    varshape = []
    vardim_sub = []
    for idim in vardim:
#        varshape.append(scale_dimdef['len'][idim])
        varshape.append(scale_dimdef['len'][idim][0])
        vardim_sub.append(None)
        for idiml in scale_dimlist_g:
            if idim in idiml:
                varshape[-1] = scale_dimdef['len_g'][idim]
                vardim_sub[-1] = idim
                break

    if all(i is None for i in vardim_sub):
        return vardim, vardata_0
    else:
        if type(vardata_0) == ma.MaskedArray:
            vardata = ma.masked_all(varshape, dtype=vardata_0.dtype)
            vardata.fill_value = vardata_0.fill_value
        else:
            vardata = np.empty(varshape, dtype=vardata_0.dtype)
        for ip in range(nproc):
            slice_obj = [slice(None)] * len(vardim)
            for i, idim in enumerate(vardim_sub):
                if idim is not None:
#                    slice_obj[i] = slice(scale_dimdef['start'][idim][ip],
#                                         scale_dimdef['start'][idim][ip] + scale_dimdef['len'][idim])
                    slice_obj[i] = slice(scale_dimdef['start'][idim][ip],
                                         scale_dimdef['start'][idim][ip] + scale_dimdef['len'][idim][ip])
            if ip == 0:
                vardata[slice_obj] = vardata_0
            else:
                vardim, vardata[slice_obj] = ncphysio.ncphys_read(rootgrps[ip], varname, dimlist=scale_dimlist, time=time, it=it)
        return vardim, vardata


def scale_write(nproc, rootgrps, scale_dimdef, varname, vardata, t=None):
    """
    Write a variable to a set of split SCALE files.
    Assume the input dimensions are consistent.

    Parameters
    ----------
    nproc : integer
        Number of split SCALE files
    rootgrps : array of netcdf4-python Dataset instance
        Array of netcdf4-python Dataset instances of the split SCALE files
    scale_dimdef : dictionary
        Summary of dimensions in the split SCALE files
    varname : string
        The variable name.
    vardata : ndarray or masked_array
        Variable data to be written to the files
    t : int or <datetime.datetime> class or None, optional
        Time to read. None for all times. Defalut: None
    """
    if t is None:
        it = 'all'
        time = 'all'
    elif type(t) is int:
        it = t
        time = None
    elif type(t) is dt.datetime:
        it = None
        time = scale_puttime(t)
    else:
        raise ValueError("The type of 't' should be either 'int' or 'datetime.datetime' or 'None'.")

    # Check if the variable exists in the NetCDF file
    if varname in rootgrps[0].variables:
        vardim_in = rootgrps[0].variables[varname].dimensions
    else:
        raise IOError("Variable '" + varname + "' does not exist.")

    vardim = []
    for i, idiml in enumerate(scale_dimlist):
        if i > 0 or (time == 'all' or it == 'all'):
            for idim in vardim_in:
                if idim in idiml:
                    vardim.append(idim)

    varshape = []
    vardim_sub = []
    for idim in vardim:
#        varshape.append(scale_dimdef['len'][idim])
        varshape.append(scale_dimdef['len'][idim][0])
        vardim_sub.append(None)
        for idiml in scale_dimlist_g:
            if idim in idiml:
                varshape[-1] = scale_dimdef['len_g'][idim]
                vardim_sub[-1] = idim
                break

    for ip in range(nproc):
        slice_obj = [slice(None)] * len(vardim)
        for i, idim in enumerate(vardim_sub):
            if idim is not None:
#                slice_obj[i] = slice(scale_dimdef['start'][idim][ip],
#                                     scale_dimdef['start'][idim][ip] + scale_dimdef['len'][idim])
                slice_obj[i] = slice(scale_dimdef['start'][idim][ip],
                                     scale_dimdef['start'][idim][ip] + scale_dimdef['len'][idim][ip])
        ncphysio.ncphys_write(rootgrps[ip], varname, vardim, vardata[slice_obj], dimlist=scale_dimlist, time=time, it=it)


class ScaleIO:
    """
    Class for split SCALE I/O

    Attributes
    ----------
    ***
    """
    def __init__(self, basename, mode='r'):
        """
        Parameters
        ----------
        basename : string
            Split SCALE file basename. Path can be included.
        mode : {'r', 'r+'}, optional
            File I/O mode
            * 'r' -- read (default)
            * 'r+' -- read/write
        """
        self.nproc, self.rootgrps, self.dimdef = scale_open(basename, mode)
        self.z = scale_read(self.nproc, self.rootgrps, self.dimdef, 'z')[1]
        if self.dimdef['len']['time'][0] is None:
            self.t = None
        else:
            time_array = scale_read(self.nproc, self.rootgrps, self.dimdef, 'time')[1]
            self.t = np.empty_like(time_array, dtype='O')
#            self.t = [scale_gettime(i) for i in time_array]
            for it in range(len(time_array)):
                self.t[it] = scale_gettime(time_array[it])
        self.z = scale_read(self.nproc, self.rootgrps, self.dimdef, 'z')[1]
        self.zh = scale_read(self.nproc, self.rootgrps, self.dimdef, 'zh')[1]
        self.lon = scale_read(self.nproc, self.rootgrps, self.dimdef, 'lon')[1]
        self.lat = scale_read(self.nproc, self.rootgrps, self.dimdef, 'lat')[1]


    def __del__(self):
        try:
            scale_close(self.rootgrps)
        except:
            pass


    def readvar(self, varname, t=None):
        """
        Read a variable from a set of split SCALE files.

        Parameters
        ----------
        varname : string
            The variable name.
        t : int or <datetime.datetime> class or None, optional
            Time to read
            * None -- all times (defalut)

        Returns
        -------
        vardata : ndarray or masked_array
            Variable data in a ndarray or masked_array (if the variable has the 
            `_FillValue` attribute).
        """
        return scale_read(self.nproc, self.rootgrps, self.dimdef, varname, t=t)[1]


    def writevar(self, varname, vardata, t=None):
        """
        Write a variable to a set of split SCALE files.
        Assume the input dimensions are consistent.

        Parameters
        ----------
        varname : string
            The variable name.
        vardata : ndarray or masked_array
            Variable data to be written to the files.
        t : int or <datetime.datetime> class or None, optional
            Time to read. None for all times. Defalut: None
        """
        scale_write(self.nproc, self.rootgrps, self.dimdef, varname, vardata, t=t)
