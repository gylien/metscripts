import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset

dimlist_default = [['time'], ['z'], ['y'], ['x']]


def ncphys_read_dimlen(rootgrp, dimlist=dimlist_default):
    """
    Read lengths of dimensions.

    Parameters
    ----------
    rootgrp : netcdf4-python Dataset instance
        The input NetCDF file.
    dimlist : array of array, optional
        List of dimensions in the NetCDF file. Default: `dimlist_default`

    Returns
    -------
    dimlen : dictionary
        Lengths of dimensions in the NetCDF file.
    """
    dimlen = {}
    for idiml in dimlist:
        for idim in idiml:
            if idim in rootgrp.dimensions:
                dimlen[idim] = len(rootgrp.dimensions[idim])
            else:
                dimlen[idim] = None
    return dimlen


def ncphys_read(rootgrp, varname, dimlist=dimlist_default, dimlen=None, time=None, it=None):
    """
    Read a variable from a NetCDF file.

    Can choose to read a single time or all times.
    The dimensions of the variable are reordered based on the 'dimlist' setting.

    Parameters
    ----------
    rootgrp : netcdf4-python Dataset instance
        The input NetCDF file.
    varname : string
        The variable name.
    dimlist : array of array, optional
        List of dimensions in the NetCDF file. Default: `dimlist_default`
    dimlen : dictionary, optional
        Lengths of dimensions in the NetCDF file, returned from 
        `ncphys_read_dims`. Default: automatically get the lengths of dimensions.
    time : number, optional
        Target time in physical time unit.
    it : int, optional
        Target time index. Defalut to 0 if both `time` and `it` are not given.

    Returns
    -------
    vardim : dictionary
        Dimensions of the return variable data.
    vardata : ndarray or masked_array
        Variable data in a ndarray or masked_array (if the variable has the 
        `_FillValue` attribute).
    """
    # If dimlen is not given, read the dimension lengths
    if dimlen is None:
        dimlen = ncphys_read_dims(rootgrp, dimlist=dimlist)

    # Check the variable dimensions and the order of the variable dimensions
    # in the input NetCDF file
    vardim = rootgrp.variables[varname].dimensions
    dimord = [None] * len(dimlist)
    dimname = [None] * len(dimlist)
    for i, idiml in enumerate(dimlist):
        for j, idim in enumerate(vardim):
            if idim in idiml:
                if dimord[i] is None:
                    dimord[i] = j
                    dimname[i] = idim
                else:
                    raise ValueError, 'Duplicated dimensions.'
    if len([i for i in dimord if i is not None]) != len(vardim):
        raise ValueError, "Variable has addtional dimensions not listed in 'dimlist'"

    # Read the variable data into a ndarray or masked_array
    if dimord[0] is None or time == 'all' or it == 'all':
        vardata = rootgrp.variables[varname][:]
    else:
        if time is not None:
            found = np.where(rootgrp.variables[dimname[0]][:] == time)[0]
            if len(found) == 0:
                raise ValueError, 'Cannot find time = ' + str(time)
            else:
                it = found[0]
        elif it is None:
            it = 0
        if dimord[0] == 0:
            vardata = rootgrp.variables[varname][it]
        else:
            vardata = np.take(rootgrp.variables[varname][:], it, axis=dimord[0])
    if '_FillValue' in rootgrp.variables[varname].ncattrs() and type(vardata) != ma.MaskedArray:
        vardata = ma.array(vardata, fill_value=rootgrp.variables[varname].getncattr('_FillValue'))

    # Reorder the axes of the return variable data if necessary
    if time == 'all' or it == 'all':
        vardim_out = [i for i in dimname if i is not None]
        transaxes = [i for i in dimord if i is not None]
    else:
        vardim_out = [i for i in dimname[1:] if i is not None]
        transaxes = []
        for iord in dimord[1:]:
            if iord is not None:
                if dimord[0] is not None and iord > dimord[0]:
                    transaxes.append(iord - 1)
                else:
                    transaxes.append(iord)
    if transaxes != range(len(transaxes)):
        vardata = np.transpose(vardata, axes=transaxes)

    return vardim_out, vardata
