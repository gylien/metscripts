import numpy as np
import numpy.ma as ma
from ncphysio import *
from netCDF4 import Dataset
import os.path
import datetime as dt

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
    Open a set of partial scale I/O files and return the definition of 
    global and subdomain dimensions.
    """
    ncfile = basename + scale_file_suffix.format(0)
    if not os.path.isfile(ncfile):
        if mode == 'r':
            raise IOError, 'File does not exist.'
        elif scale_dimdef is None:
            raise IOError, 'File does not exist.'
        else:
#            scale_create_new(basename, scale_dimdef)
            raise IOError, 'Scale_create has not been supported yet...'

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
                    sub_idx[idim] = range(len(rootgrps[ip].dimensions[idim]))
                    sub_var[idim] = rootgrps[ip].variables[idim][:]
                else:
                    sub_ip[idim] += [ip] * len(rootgrps[ip].dimensions[idim])
                    sub_idx[idim] += range(len(rootgrps[ip].dimensions[idim]))
                    sub_var[idim] = np.append(sub_var[idim], rootgrps[ip].variables[idim][:])
        ip += 1
    nproc = ip
    dimlen = ncphys_read_dimlen(rootgrps[0], dimlist=scale_dimlist)

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
                    raise ValueError, 'Subdomains are not consistent.'

    scale_dimdef = {'len': dimlen, 'len_g': dimlen_g, 'coor_g': dimcoor_g, 'start': dimstart}
    return nproc, rootgrps, scale_dimdef


def scale_close(rootgrps):
    """
    Close a set of partial scale I/O files.
    """
    for irg in rootgrps:
        irg.close()


def scale_read(nproc, rootgrps, scale_dimdef, varname, time=None, it=None):
    """
    Read a variable from a set of partial scale I/O files.

    Can choose to read a single time or all times.
    """
    vardim, vardata_0 = ncphys_read(rootgrps[0], varname, dimlist=scale_dimlist, time=time, it=it)
    varshape = []
    vardim_sub = []
    for idim in vardim:
        varshape.append(scale_dimdef['len'][idim])
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
        for ip in xrange(nproc):
            slice_obj = [slice(None)] * len(vardim)
            for i, idim in enumerate(vardim_sub):
                if idim is not None:
                    slice_obj[i] = slice(scale_dimdef['start'][idim][ip],
                                         scale_dimdef['start'][idim][ip] + scale_dimdef['len'][idim])
            if ip == 0:
                vardata[slice_obj] = vardata_0
            else:
                vardim, vardata[slice_obj] = ncphys_read(rootgrps[ip], varname, dimlist=scale_dimlist, time=time, it=it)
        return vardim, vardata


def scale_write(nproc, rootgrps, scale_dimdef, varname, vardim, vardata, time=None, it=None):
    """
    Write a variable to a set of partial scale I/O files.

    Can choose to write a single time or all times.
    """
    varshape = []
    vardim_sub = []
    for idim in vardim:
        varshape.append(scale_dimdef['len'][idim])
        vardim_sub.append(None)
        for idiml in scale_dimlist_g:
            if idim in idiml:
                varshape[-1] = scale_dimdef['len_g'][idim]
                vardim_sub[-1] = idim
                break

    for ip in xrange(nproc):
        slice_obj = [slice(None)] * len(vardim)
        for i, idim in enumerate(vardim_sub):
            if idim is not None:
                slice_obj[i] = slice(scale_dimdef['start'][idim][ip],
                                     scale_dimdef['start'][idim][ip] + scale_dimdef['len'][idim])
        ncphys_write(rootgrps[ip], varname, vardim, vardata[slice_obj], dimlist=scale_dimlist, time=time, it=it)


def scale_gettime(scale_time):
    """
    Convert the python datetime to the scale model time.
    """
    return scale_time_0 + dt.timedelta(seconds=scale_time)


def scale_puttime(datetime):
    """
    Convert the scale model time to the python datetime.
    """
    return (datetime - scale_time_0).total_seconds()
