import numpy as np
import numpy.ma as ma


#__all__ = ['scale_dimlist', 'scale_dimlist_g', 'scale_file_suffix', 'scale_time_0',
#           'scale_open', 'scale_close', 'scale_gettime', 'scale_puttime',
#           'scale_read', 'scale_write', 'ScaleIO']


def calc_height(sio, topo=None):
    """
    Calculate the 3-D height

    Parameters
    ----------
    sio : <scale.io.ScaleIO> class
        Split SCALE I/O class
    topo : 2-D ndarray, optional
        Surface height (m)

    Returns
    -------
    height : 3-D ndarray
        Height in full levels (m)
    height_h : 3-D ndarray
        Height in half levels (m)
    """
    if topo is None:
        topo0 = sio.readvar('TOPO')
    else:
        topo0 = topo

    height = np.zeros((sio.dimdef['len']['z'][0], topo0.shape[0], topo0.shape[1]), dtype=topo0.dtype)
    height_h = np.zeros((sio.dimdef['len']['zh'][0], topo0.shape[0], topo0.shape[1]), dtype=topo0.dtype)

    for k in range(len(sio.z)):
        height[k,:,:] = topo0 + (sio.zh[-1] - topo0) / sio.zh[-1] * sio.z[k]
    for k in range(len(sio.zh)):
        height_h[k,:,:] = topo0 + (sio.zh[-1] - topo0) / sio.zh[-1] * sio.zh[k]

    return height, height_h


def interp_z(sio, var, height=None):
    """
    """
    if height is None:
        height0 = sio.readvar('height')
    else:
        height0 = height
    height0[0,:,:] -= 1.e-2

    varout = ma.masked_all((sio.dimdef['len']['z'][0], var.shape[1], var.shape[2]), dtype=var.dtype)
    if type(var) == ma.MaskedArray:
        varout.fill_value = var.fill_value

    for j in range(var.shape[1]):
        for i in range(var.shape[2]):
            varout[:,j,i] = np.interp(sio.z, height0[:,j,i], var[:,j,i], left=-9.9999e+30, right=-9.9999e+30)
    varout.mask = [varout == -9.9999e+30]

    return varout



