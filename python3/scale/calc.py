import numpy as np
import numpy.ma as ma


__all__ = ['calc_height', 'interp_z', 'interp_p', 'calc_destagger', 'calc_pt', 'calc_uvw', 'calc_ref']


def calc_height(sio, topo=None):
    """
    Calculate the 3-D height

    Parameters
    ----------
    sio : <scale.io.ScaleIO> class
        Split SCALE I/O class
    topo : 2-D ndarray, optional
        Surface height (m). Read from files if not given

    Returns
    -------
    height : 3-D ndarray
        Height in full levels (m)
    height_h : 3-D ndarray
        Height in half levels (m)
    """
    if topo is None:
        topo0 = sio.readvar('TOPO')[2:-2,2:-2]
    else:
        topo0 = topo

    height = np.zeros((sio.dimdef['len']['z'][0], topo0.shape[0], topo0.shape[1]), dtype=topo0.dtype)
    height_h = np.zeros((sio.dimdef['len']['zh'][0], topo0.shape[0], topo0.shape[1]), dtype=topo0.dtype)

    for k in range(len(sio.z)):
        height[k,:,:] = topo0 + (sio.zh[-1] - topo0) / sio.zh[-1] * sio.z[k]
    for k in range(len(sio.zh)):
        height_h[k,:,:] = topo0 + (sio.zh[-1] - topo0) / sio.zh[-1] * sio.zh[k]

    return height, height_h


def interp_z(sio, var, height=None, t=None):
    """
    Interpolate a 3-D variable to constant z levels

    Parameters
    ----------
    sio : <scale.io.ScaleIO> class
        Split SCALE I/O class
    var : 3-D ndarray
        Input variable
    height : 3-D ndarray, optional
        Height in current grids (m). Read from files if not given
    t : int or <datetime.datetime> class or None, optional
        Time to read
        * None -- all times (defalut)

    Returns
    -------
    varout : 3-D ndarray
        Interpolated variable
    """
    if height is None:
        height0 = sio.readvar('height', t=t)[:,2:-2,2:-2]
    else:
        height0 = height
    height0[0,:,:] -= 1.e-2

    varshape = list(var.shape)
    varshape[0] = sio.dimdef['len']['z'][0]
    varout = ma.masked_all(varshape, dtype=var.dtype)
    if type(var) == ma.MaskedArray:
        varout.fill_value = var.fill_value

    for j in range(varshape[1]):
        for i in range(varshape[2]):
            varout[:,j,i] = np.interp(sio.z, height0[:,j,i], var[:,j,i], left=-9.99e+33, right=-9.99e+33)
    varout.mask[varout == -9.99e+33] = True

    return varout


def interp_p(sio, var, plevels, p=None, t=None):
    """
    Interpolate a 3-D variable to constant p levels

    Parameters
    ----------
    sio : <scale.io.ScaleIO> class
        Split SCALE I/O class
    var : 3-D ndarray
        Input variable
    plevels : 1-D ndarray
        Targeted pressure levels (Pa)
    p : 3-D ndarray, optional
        Pressure in current grids (Pa). Read from files if not given
    t : int or <datetime.datetime> class or None, optional
        Time to read
        * None -- all times (defalut)

    Returns
    -------
    varout : 3-D ndarray
        Interpolated variable
    """
    if p is None:
        p0, = calc_pt(sio, tout=False, thetaout=False, t=t)
    else:
        p0 = p
    varshape = list(var.shape)
    varshape[0] = len(plevels)
    varout = ma.masked_all(varshape, dtype=var.dtype)
    if type(var) == ma.MaskedArray:
        varout.fill_value = var.fill_value

    for j in range(varshape[1]):
        for i in range(varshape[2]):
            varout[:,j,i] = np.interp(np.log(plevels), np.log(p0[::-1,j,i]), var[::-1,j,i], left=-9.99e+33, right=-9.99e+33)
    varout.mask[varout == -9.99e+33] = True

    return varout


def calc_destagger(var, axis=0, first_grd=False):
    """
    Calculate full-level values from the half-level (staggered-grid) values

    Parameters
    ----------
    var : ndarray
        Input variable
    axis : int, optional
        Staggered axis. Default: 0
    first_grd : bool, optional
        * True -- Addtional first-row grids are provided for interpolation
        * False -- No additional first-row grid (default)

    Returns
    -------
    varout : 3-D ndarray
        Destaggered variable
    """
    if axis < 0 or axis >= var.ndim:
        raise ValueError("'axis' is invalid. It must be within [0, var.ndim).")
    varshape = list(var.shape)
    if first_grd:
        varshape[axis] -= 1
    if type(var) == ma.MaskedArray:
        varout = ma.masked_all(varshape, dtype=var.dtype)
        varout.fill_value = var.fill_value
    else:
        varout = np.empty(varshape, dtype=var.dtype)

    slice_obj = [slice(None)] * var.ndim
    if first_grd:
        for i in range(varshape[axis]):
            slice_obj[axis] = i
            varout[slice_obj] = 0.5 * (var.take(i, axis=axis) + var.take(i+1, axis=axis))
    else:
        slice_obj[axis] = 0
        varout[slice_obj] = var.take(0, axis=axis)
        for i in range(1, varshape[axis]):
            slice_obj[axis] = i
            varout[slice_obj] = 0.5 * (var.take(i-1, axis=axis) + var.take(i, axis=axis))

    return varout


def calc_pt(sio, rho=None, rhot=None, qv=None, qhydro=None, tout=True, thetaout=False, t=None):
    """
    Calculate 3-D pressure, temperature, potential temperature

    Parameters
    ----------
    sio : <scale.io.ScaleIO> class
        Split SCALE I/O class
    rho : 3-D ndarray, optional
        Density (kg/m3). Read from files if not given
    rhot : 3-D ndarray, optional
        rho * theta (kg/m3*K). Read from files if not given
    qv : 3-D ndarray, optional
        Water vapor mixing ratio (kg/kg). Read from files if not given
    qhydro : 3-D ndarray, optional
        Mixing ratio of all hydrometers (kg/kg). Read from files if not given
    tout : bool, optional
        Whether output temperature? Default: True
    thetaout : bool, optional
        Whether output potential temperature? Default: False
    t : int or <datetime.datetime> class or None, optional
        Time to read
        * None -- all times (defalut)

    Returns
    -------
    p : 3-D ndarray
        Pressure (Pa)
    t : 3-D ndarray, optional
        Temperature (K)
    theta : 3-D ndarray, optional
        Potential temperature (K)
    """
    if rho is None:
        rho0 = sio.readvar('DENS', t=t)[:,2:-2,2:-2]
    else:
        rho0 = rho
    if rhot is None:
        rhot0 = sio.readvar('RHOT', t=t)[:,2:-2,2:-2]
    else:
        rhot0 = rhot
    if qv is None:
        qv0 = sio.readvar('QV', t=t)[:,2:-2,2:-2]
    else:
        qv0 = qv
    if qhydro is None:
        qhydro0 = sio.readvar('QC', t=t)[:,2:-2,2:-2] \
                + sio.readvar('QR', t=t)[:,2:-2,2:-2] \
                + sio.readvar('QI', t=t)[:,2:-2,2:-2] \
                + sio.readvar('QS', t=t)[:,2:-2,2:-2] \
                + sio.readvar('QG', t=t)[:,2:-2,2:-2]
    else:
        qhydro0 = qhydro

    Rdry = 287.04
    Rvap = 461.46
    CVdry = 717.60
    PRE00 = 100000.0

    Rtot = Rdry * (1. - qv0 - qhydro0) + Rvap * qv0
    CPovCV = ( CVdry + Rtot ) / CVdry

    p = PRE00 * np.power(rhot0 * Rtot / PRE00, CPovCV)
    res = [p]
    if tout:
        t = p / (rho0 * Rtot)
        res.append(t)
    if thetaout:
        theta = rhot0 / rho0
        res.append(theta)

    return res



def calc_uvw(sio, rho=None, momx=None, momy=None, momz=None, destagger=True, first_grd=True, t=None):
    """
    Calculate 3-D u, v, w winds

    Parameters
    ----------
    sio : <scale.io.ScaleIO> class
        Split SCALE I/O class
    rho : 3-D ndarray, optional
        Density (kg/m3). Read from files if not given
    momx : 3-D ndarray, optional
        x-momentum (kg/m2/s). Read from files if not given
    momy : 3-D ndarray, optional
        y-momentum (kg/m2/s). Read from files if not given
    momz : 3-D ndarray, optional
        z-momentum (kg/m2/s). Read from files if not given
    destaggered : bool, optional
        * True -- Destagger momx, momy, momz before calculation (default)
        * False -- Do not need to destagger momx, momy, momz
    first_grd : bool, optional
        * True -- Addtional first-row grids are provided for interpolation (default)
        * False -- No additional first-row grid
    t : int or <datetime.datetime> class or None, optional
        Time to read
        * None -- all times (defalut)

    Returns
    -------
    u : 3-D ndarray
        u-wind (m/s)
    v : 3-D ndarray
        v-wind (m/s)
    w : 3-D ndarray
        w-wind (m/s)
    """
    if rho is None:
        rho0 = sio.readvar('DENS', t=t)[:,2:-2,2:-2]
    else:
        rho0 = rho
    if momx is None:
        if first_grd:
            momx0 = sio.readvar('MOMX', t=t)[:,2:-2,1:-2]
        else:
            momx0 = sio.readvar('MOMX', t=t)[:,2:-2,2:-2]
    else:
        momx0 = momx
    if momy is None:
        if first_grd:
            momy0 = sio.readvar('MOMY', t=t)[:,1:-2,2:-2]
        else:
            momy0 = sio.readvar('MOMY', t=t)[:,2:-2,2:-2]
    else:
        momy0 = momy
    if momz is None:
        momz0 = sio.readvar('MOMZ', t=t)[:,2:-2,2:-2]
    else:
        momz0 = momz

    if destagger:
        momx0 = calc_destagger(momx0, axis=2, first_grd=first_grd)
        momy0 = calc_destagger(momy0, axis=1, first_grd=first_grd)
        momz0 = calc_destagger(momz0, axis=0, first_grd=False)
    u = momx0 / rho0
    v = momy0 / rho0
    w = momz0 / rho0
    w[0,:,:] = 0.

    return u, v, w


def calc_ref(sio, rho=None, qr=None, qs=None, qg=None, t=None):
    """
    Calculate radar reflectivity

    Parameters
    ----------
    sio : <scale.io.ScaleIO> class
        Split SCALE I/O class
    rho : 3-D ndarray, optional
        Density (kg/m3). Read from files if not given
    qr : 3-D ndarray, optional
        Rain water mixing ratio (kg/kg). Read from files if not given
    qs : 3-D ndarray, optional
        Snow mixing ratio (kg/kg). Read from files if not given
    qg : 3-D ndarray, optional
        Graupel mixing ratio (kg/kg). Read from files if not given
    t : int or <datetime.datetime> class or None, optional
        Time to read
        * None -- all times (defalut)

    Returns
    -------
    dbz : 3-D ndarray
        Radar reflectivity (dBZ)
    max_dbz : 2-D ndarray
        Maximum radar reflectivity in vertical (dBZ)
    """
    if rho is None:
        rho0 = sio.readvar('DENS', t=t)[:,2:-2,2:-2]
    else:
        rho0 = rho
    if qr is None:
        qr0 = sio.readvar('QR', t=t)[:,2:-2,2:-2]
    else:
        qr0 = qr
    if qs is None:
        qs0 = sio.readvar('QS', t=t)[:,2:-2,2:-2]
    else:
        qs0 = qs
    if qg is None:
        qg0 = sio.readvar('QG', t=t)[:,2:-2,2:-2]
    else:
        qg0 = qg

    qr0[qr0 < 0.] = 0.
    qs0[qs0 < 0.] = 0.
    qg0[qg0 < 0.] = 0.
    ref = 2.53e4 * (rho0 * qr0 * 1.0e3) ** 1.84 \
        + 3.48e3 * (rho0 * qs0 * 1.0e3) ** 1.66 \
        + 8.18e4 * (rho0 * qg0 * 1.0e3) ** 1.50
    dbz = 10. * np.log10(ref)
    max_dbz = ma.max(dbz, axis=0)

    min_dbz = -20.
    dbz[dbz < min_dbz] = min_dbz
    max_dbz[max_dbz < min_dbz] = min_dbz

    return dbz, max_dbz
