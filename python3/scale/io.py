import numpy as np
import numpy.ma as ma
import ncphysio
import os.path
import datetime as dt
from netCDF4 import Dataset


__all__ = ['scale_dimlist', 'scale_dimlist_g', 'scale_file_suffix',
           'scale_open', 'scale_close', 'scale_gettime', 'scale_puttime',
           'scale_read', 'scale_write', 'ScaleIO']


show_cache_use = True


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

###
###        print(ncfile)
###

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


def scale_gettime(scale_time, year):
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
    return dt.datetime(year, 1, 1, 0, 0, 0) + dt.timedelta(seconds=scale_time)


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
    return (time - dt.datetime(time.year, 1, 1, 0, 0, 0)).total_seconds()


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

###
###    print(varshape)
###


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

###
###            print(varname, ip)
###


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
    def __init__(self, basename, mode='r', year=None, bufsize=0, cache=False):
        """
        Parameters
        ----------
        basename : string
            Split SCALE file basename. Path can be included.
        mode : {'r', 'r+'}, optional
            File I/O mode
            * 'r' -- read (default)
            * 'r+' -- read/write
        year : integer
            Year of the data.
            * None -- the current year (default)
        bufsize : int
            Unused grid numbers near the lateral boundary (default: 0)
        cache : bool
            Option to cache the data read into the memory.
            * False -- do not cache (default)
            * True -- cache
        """
        if year is None:
            self.year = dt.datetime.now().year
        self.nproc, self.rootgrps, self.dimdef = scale_open(basename, mode)
        self.z = scale_read(self.nproc, self.rootgrps, self.dimdef, 'z')[1]
        if self.dimdef['len']['time'][0] is None:
            self.t = None
        else:
            time_array = scale_read(self.nproc, self.rootgrps, self.dimdef, 'time')[1]
            self.t = np.empty_like(time_array, dtype='O')
#            self.t = [scale_gettime(i, self.year) for i in time_array]
            for it in range(len(time_array)):
                self.t[it] = scale_gettime(time_array[it], self.year)
        self.z = scale_read(self.nproc, self.rootgrps, self.dimdef, 'z')[1]
        self.zh = scale_read(self.nproc, self.rootgrps, self.dimdef, 'zh')[1]
        self.lon = scale_read(self.nproc, self.rootgrps, self.dimdef, 'lon')[1]
        self.lat = scale_read(self.nproc, self.rootgrps, self.dimdef, 'lat')[1]
        assert bufsize >= 0, "'bufsize' should be greater than or equal to 0."
        self.bufsize = bufsize
        if cache:
            self.cache = {}


    def __del__(self):
        try:
            scale_close(self.rootgrps)
        except:
            pass


    def freecache(self):
        if hasattr(self, 'cache'):
            del self.cache
            self.cache = {}
        else:
            import warnings
            warnings.warn('Cache is not enabled.')


    def readvar(self, varname, t=None, bufsize=None):
        """
        Read a variable from a set of split SCALE files.

        Parameters
        ----------
        varname : string
            The variable name.
        t : int or <datetime.datetime> class or None, optional
            Time to read
            * None -- all times (defalut)
        bufsize : int
            Unused grid numbers near the lateral boundary
            * None -- use the bufsize given at the object initialization
                      (default)

        Returns
        -------
        vardata : ndarray or masked_array
            Variable data in a ndarray or masked_array (if the variable has the 
            `_FillValue` attribute).
        """
        if hasattr(self, 'cache'):
            if self.t is None:
                tkey = 0
            elif t is None:
                tkey = 'all'
            elif type(t) is int:
                tkey = t
            elif type(t) is dt.datetime:
                try:
                    tkey = self.t.index(t)
                except ValueError:
                    raise ValueError("Cannot find 't' = " + str(t))
            else:
                raise ValueError("The type of 't' should be either 'int' or 'datetime.datetime' or 'None'.")

            if varname in self.cache:
                if tkey in self.cache[varname]:
                    if show_cache_use:
                        if tkey == 'all':
                            print('*** read from cache: {:s} ***'.format(varname))
                        else:
                            print('*** read from cache: {:s}[t={:s}] ***'.format(varname, str(tkey)))
                    res = self.cache[varname][tkey]
                elif 'all' in self.cache[varname] and tkey != 'all':
                    if show_cache_use:
                        print('*** read from cache: {:s}[t={:s}] ***'.format(varname, str(tkey)))
                    res = self.cache[varname]['all'][tkey]
                else:
                    if tkey == 'all':
                        del self.cache[varname]
                        self.cache[varname] = {}
                    self.cache[varname][tkey] = scale_read(self.nproc, self.rootgrps, self.dimdef, varname, t=t)[1]
                    res = self.cache[varname][tkey]
            else:
                self.cache[varname] = {}
                self.cache[varname][tkey] = scale_read(self.nproc, self.rootgrps, self.dimdef, varname, t=t)[1]
                res = self.cache[varname][tkey]
        else:
            res = scale_read(self.nproc, self.rootgrps, self.dimdef, varname, t=t)[1]

        if bufsize is None:
            bufsize = self.bufsize
        assert bufsize >= 0, "'bufsize' should be greater than or equal to 0."
        if bufsize == 0 or len(res.shape) < 2:
            return res
        else:
            return res[[slice(None)] * (len(res.shape)-2) + [slice(bufsize, -bufsize), slice(bufsize, -bufsize)]]


    def writevar(self, varname, vardata, t=None, bufsize=None):
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
        bufsize : int
            Unused grid numbers near the lateral boundary
            * None -- use the bufsize given at the object initialization
                      (default)
        """
        if bufsize is None:
            bufsize = self.bufsize
        assert bufsize >= 0, "'bufsize' should be greater than or equal to 0."
        if bufsize == 0 or len(res.shape) < 2:
            scale_write(self.nproc, self.rootgrps, self.dimdef, varname, vardata, t=t)
        else:
            tmpdata = scale_read(self.nproc, self.rootgrps, self.dimdef, varname, t=t)[1]
            tmpdata[[slice(None)] * (len(tmpdata.shape)-2) + [slice(bufsize, -bufsize), slice(bufsize, -bufsize)]] = vardata
            scale_write(self.nproc, self.rootgrps, self.dimdef, varname, tmpdata, t=t)
