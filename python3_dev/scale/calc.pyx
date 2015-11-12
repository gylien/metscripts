print('hello from cython')

import numpy as np
import numpy.ma as ma
import numpy.testing as npt
from mpl_toolkits.basemap import Basemap

cimport numpy as np
cimport cython
ctypedef fused f_real:
    cython.float
    cython.double


           
__all__ = ['set_bmap', 'calc_rotate_winds',
           'calc_height', 'interp_z', 'interp_p', 'calc_destagger', 'calc_destagger_uvw', 'calc_qhydro', 'calc_pt',
           'calc_ref', 'extrap_z_t0', 'extrap_z_pt', 'extrap_p_zt', 'calc_slp', 'calc_rhosfc_psfc']


missingv = 1.e-33
rsphere = 6.37122e6


def set_bmap(sio, proj):
    """
    Set map projection

    XXXXXX
    """
    llcrnrlon = sio.lon[0, 0]
    llcrnrlat = sio.lat[0, 0]
    urcrnrlon = sio.lon[-1, -1]
    urcrnrlat = sio.lat[-1, -1]
    if proj['type'] == 'LC':
        bmap = Basemap(projection='lcc', lat_1=proj['LC_lat1'], lat_2=proj['LC_lat2'], lon_0=proj['basepoint_lon'],
                       llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,
                       rsphere=rsphere)
    else:
        raise ValueError('[Error] Unsupport map projection.')

    # verify projection setting
    x, y = bmap(sio.lon, sio.lat)
    x += sio.dimdef['coor_g']['x'][0]
    y += sio.dimdef['coor_g']['y'][0]
    for iy in range(sio.dimdef['len_g']['y']):
        npt.assert_almost_equal(x[iy,:], sio.dimdef['coor_g']['x'], decimal=6, err_msg='[Error] Incorrect projection settings.')
    for ix in range(sio.dimdef['len_g']['x']):
        npt.assert_almost_equal(y[:,ix], sio.dimdef['coor_g']['y'], decimal=6, err_msg='[Error] Incorrect projection settings.')

    return bmap


def calc_rotate_winds(sio, bmap, u=None, v=None, t=None, dryrun=False):
    """
    Calculate the rotation of u, v winds

    Parameters
    ----------
    sio : <scale.io.ScaleIO> class
        Split SCALE I/O class
    bmap : <mpl_toolkits.basemap.Basemap> class
        Basemap class
    u : 3-D ndarray, optional
        Grid wind along i direction (m/s). Read from files if not given
    v : 3-D ndarray, optional
        Grid wind along j direction (m/s). Read from files if not given
    t : int or <datetime.datetime> class or None, optional
        Time to read
        * None -- all times (defalut)
    dryrun : bool, optional
        * True -- dry run mode, only reading necessary data using 'sio' class
        * False -- do real computation (default)

    Returns
    -------
    u_rot : 3-D ndarray
        Rotated u-wind (m/s)
    v_rot : 3-D ndarray
        Rotated v-wind (m/s)
    """
    u_ = (sio.readvar('U') if u is None else u)
    v_ = (sio.readvar('V') if v is None else v)
    if dryrun:
        return None, None

    if sio.bufsize == 0:
        lon = np.copy(sio.lon)
        lat = np.copy(sio.lat)
    else:
        lon = np.copy(sio.lon[sio.bufsize:-sio.bufsize, sio.bufsize:-sio.bufsize])
        lat = np.copy(sio.lat[sio.bufsize:-sio.bufsize, sio.bufsize:-sio.bufsize])
    ny, nx = lon.shape

    if len(u_.shape) == 3:
        nz = u_.shape[0]
        lon = np.repeat(lon[np.newaxis,:,:], nz, axis=0)
        lat = np.repeat(lat[np.newaxis,:,:], nz, axis=0)

    u_rot, v_rot = bmap.rotate_vector(u_, v_, lon, lat)

    mag_square = u_ * u_ + v_ * v_
    tmpcos = (u_rot * u_ + v_rot * v_) / mag_square
    tmpsin = (u_rot * v_ - v_rot * u_) / mag_square
    u_rot = (tmpcos * u_ - tmpsin * v_).astype(u_.dtype)
    v_rot = (tmpsin * u_ + tmpcos * v_).astype(u_.dtype)

    return u_rot, v_rot


def calc_height(sio, topo=None, dryrun=False):
    """
    Calculate the 3-D height

    Parameters
    ----------
    sio : <scale.io.ScaleIO> class
        Split SCALE I/O class
    topo : 2-D ndarray, optional
        Surface height (m). Read from files if not given
    dryrun : bool, optional
        * True -- dry run mode, only reading necessary data using 'sio' class
        * False -- do real computation (default)

    Returns
    -------
    height : 3-D ndarray
        Height in full levels (m)
    height_h : 3-D ndarray
        Height in half levels (m)
    """
    topo_ = (sio.readvar('TOPO') if topo is None else topo)
    if dryrun:
        return None, None

    height = np.zeros((sio.dimdef['len']['z'][0], topo_.shape[0], topo_.shape[1]), dtype=topo_.dtype)
    height_h = np.zeros((sio.dimdef['len']['zh'][0], topo_.shape[0], topo_.shape[1]), dtype=topo_.dtype)

    for k in range(len(sio.z)):
        height[k,:,:] = topo_ + (sio.zh[-1] - topo_) / sio.zh[-1] * sio.z[k]
    for k in range(len(sio.zh)):
        height_h[k,:,:] = topo_ + (sio.zh[-1] - topo_) / sio.zh[-1] * sio.zh[k]

    return height, height_h


def interp_z(sio, np.ndarray[f_real, ndim=3] var, height=None, t=None, extrap=False, dryrun=False):
#def interp_z(sio, var, height=None, t=None, extrap=False):
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
    extrap : bool, optional
        Extrapolate low values?
        * True -- do extrapolation
        * False -- do not do extrapolation (default)
    dryrun : bool, optional
        * True -- dry run mode, only reading necessary data using 'sio' class
        * False -- do real computation (default)

    Returns
    -------
    varout : 3-D ndarray
        Interpolated variable
    """
    cdef np.ndarray[f_real, ndim=3] height_ = (sio.readvar('height', t=t) if height is None else height)
    if dryrun:
        return None

    height_[0,:,:] -= 1.e-2

    cdef int nk = var.shape[0]
    cdef int nj = var.shape[1]
    cdef int ni = var.shape[2]
    cdef int nko = sio.dimdef['len']['z'][0]
#    varshape = list(var.shape)
#    varshape[0] = sio.dimdef['len']['z'][0]
    cdef np.ndarray[double, ndim=3] varout = np.empty((nko,nj,ni), dtype='f8')

    cdef np.ndarray[double, ndim=1] zlevels = sio.z

    cdef int i, j, k, k1, ko
    for j in range(nj):
        for i in range(ni):

            k1 = 0
            for ko in range(nko):
                for k in range(k1, nk+1):
                    if k == nk:
                        break
                    elif height_[k,j,i] > zlevels[ko]: # assume an ascending order of zlevels
                        break
                k1 = k
                if k1 == 0:
                    if extrap:
                        varout[ko,j,i] = var[0,j,i]
                    else:
                        varout[ko,j,i] = missingv
                elif k1 == nk:
                    varout[ko,j,i] = missingv
                else:
                    varout[ko,j,i] = var[k1-1,j,i] + \
                                     (var[k1,j,i] - var[k1-1,j,i]) * (zlevels[ko] - height_[k1-1,j,i]) / (height_[k1,j,i] - height_[k1-1,j,i])
#            if extrap:
#                varout[:,j,i] = np.interp(sio.z, np.copy(height_[:,j,i]), np.copy(var[:,j,i]), right=missingv)
#            else:
#                varout[:,j,i] = np.interp(sio.z, np.copy(height_[:,j,i]), np.copy(var[:,j,i]), left=missingv, right=missingv)

    varout_ma = ma.masked_where(varout == missingv, varout)
    if type(var) == ma.MaskedArray:
        varout_ma.fill_value = var.fill_value

    return varout_ma


def interp_p(sio, np.ndarray[f_real, ndim=3] var, plevels, p=None, t=None, bint extrap=False, dryrun=False):
#def interp_p(sio, var, plevels, p=None, t=None, extrap=False, threads=1):
    """
    Interpolate a 3-D variable to constant p levels

    Parameters
    ----------
    sio : <scale.io.ScaleIO> class
        Split SCALE I/O class
    var : 3-D ndarray
        Input variable
    plevels : array_like, assuming a descending order (from bottom to top)
        Targeted pressure levels (Pa)
    p : 3-D ndarray, optional
        Pressure in current grids (Pa). Read from files if not given
    t : int or <datetime.datetime> class or None, optional
        Time to read
        * None -- all times (defalut)
    extrap : bool, optional
        Extrapolate low values?
        * True -- do extrapolation
        * False -- do not do extrapolation (default)
    dryrun : bool, optional
        * True -- dry run mode, only reading necessary data using 'sio' class
        * False -- do real computation (default)

    Returns
    -------
    varout : 3-D ndarray
        Interpolated variable
    """
    cdef np.ndarray[f_real, ndim=3] p_ = (calc_pt(sio, tout=False, thetaout=False, t=t, dryrun=dryrun)[0] if p is None else p)
    if dryrun:
        return None

    cdef int nk = var.shape[0]
    cdef int nj = var.shape[1]
    cdef int ni = var.shape[2]
    cdef int nko = len(plevels)
#    varshape = list(var.shape)
#    varshape[0] = len(plevels)
    cdef np.ndarray[double, ndim=3] varout = np.empty((nko,nj,ni), dtype='f8')

    cdef np.ndarray[f_real, ndim=3] log_p = np.log(p_)
    cdef np.ndarray[double, ndim=1] log_plevels = np.log(plevels)
#    cdef np.ndarray[double, ndim=3] log_p_inv = np.log(np.transpose(p_, (1, 2, 0))[:,:,::-1])
#    cdef np.ndarray[double, ndim=3] var_inv = np.copy(np.transpose(var, (1, 2, 0))[:,:,::-1])

    cdef int i, j, k, k1, ko
    for j in range(nj):
        for i in range(ni):

            k1 = 0
            for ko in range(nko):
                for k in range(k1, nk+1):
                    if k == nk:
                        break
                    elif log_p[k,j,i] < log_plevels[ko]: # assume a descending order of plevels
                        break
                k1 = k
                if k1 == 0:
                    if extrap:
                        varout[ko,j,i] = var[0,j,i]
                    else:
                        varout[ko,j,i] = missingv
                elif k1 == nk:
                    varout[ko,j,i] = missingv
                else:
                    varout[ko,j,i] = var[k1-1,j,i] + \
                                     (var[k1,j,i] - var[k1-1,j,i]) * (log_plevels[ko] - log_p[k1-1,j,i]) / (log_p[k1,j,i] - log_p[k1-1,j,i])
#            if extrap:
#                varout[:,j,i] = np.interp(log_plevels, log_p_inv[j,i,:], var_inv[j,i,:], left=missingv)
#            else:
#                varout[:,j,i] = np.interp(log_plevels, log_p_inv[j,i,:], var_inv[j,i,:], left=missingv, right=missingv)

    varout_ma = ma.masked_where(varout == missingv, varout)
    if type(var) == ma.MaskedArray:
        varout_ma.fill_value = var.fill_value

    return varout_ma


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
    varout : ndarray
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

    slice_obj_1 = [slice(None)] * var.ndim
    slice_obj_2 = [slice(None)] * var.ndim
    if first_grd:
        slice_obj_1[axis] = slice(0, varshape[axis])
        slice_obj_2[axis] = slice(1, varshape[axis]+1)
        varout = 0.5 * (var[slice_obj_1] + var[slice_obj_2])
    else:
        slice_obj_1[axis] = 0
        varout[slice_obj_1] = var[slice_obj_1]
        slice_obj_1[axis] = slice(0, varshape[axis]-1)
        slice_obj_2[axis] = slice(1, varshape[axis])
        varout[slice_obj_2] = 0.5 * (var[slice_obj_1] + var[slice_obj_2])

    return varout


def calc_destagger_uvw(sio, rho=None, momx=None, momy=None, momz=None, destagger=True, first_grd=True, t=None, dryrun=False):
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
    dryrun : bool, optional
        * True -- dry run mode, only reading necessary data using 'sio' class
        * False -- do real computation (default)

    Returns
    -------
    u : 3-D ndarray
        Destaggered u-wind (m/s)
    v : 3-D ndarray
        Destaggered v-wind (m/s)
    w : 3-D ndarray
        Destaggered w-wind (m/s)
    momx : 3-D ndarray
        Destaggered x-momentum (kg/m2/s)
    momy : 3-D ndarray
        Destaggered y-momentum (kg/m2/s)
    momz : 3-D ndarray
        Destaggered z-momentum (kg/m2/s)
    """
    rho_ = (sio.readvar('DENS', t=t) if rho is None else rho)
    momz_ = (sio.readvar('MOMZ', t=t) if momz is None else momz)
    if momx is None:
        if first_grd:
            momx_ = sio.readvar('MOMX', t=t, bufsize=0)[:,sio.bufsize:-sio.bufsize,sio.bufsize-1:-sio.bufsize]
        else:
            momx_ = sio.readvar('MOMX', t=t)
    else:
        momx_ = momx
    if momy is None:
        if first_grd:
            momy_ = sio.readvar('MOMY', t=t, bufsize=0)[:,sio.bufsize-1:-sio.bufsize,sio.bufsize:-sio.bufsize]
        else:
            momy_ = sio.readvar('MOMY', t=t)
    else:
        momy_ = momy
    if dryrun:
        return None, None, None, None, None, None

    if destagger:
#        print(' --- destagger momx')
        momx_ = calc_destagger(momx_, axis=2, first_grd=first_grd)
#        print(' --- destagger momy')
        momy_ = calc_destagger(momy_, axis=1, first_grd=first_grd)
#        print(' --- destagger momz')
        momz_ = calc_destagger(momz_, axis=0, first_grd=False)
    momz_[0,:,:] = 0.

#    print(2)

    u = momx_ / rho_
    v = momy_ / rho_
    w = momz_ / rho_

#    print(3)

    return u, v, w, momx_, momy_, momz_


def calc_qhydro(sio, qc=None, qr=None, qi=None, qs=None, qg=None, t=None, dryrun=False):
    """
    Calculate 3-D mixing ratio of all hydrometers

    Parameters
    ----------
    sio : <scale.io.ScaleIO> class
        Split SCALE I/O class
    qc : 3-D ndarray, optional
        Cloud water mixing ratio (kg/kg). Read from files if not given
    qr : 3-D ndarray, optional
        Rain mixing ratio (kg/kg). Read from files if not given
    qi : 3-D ndarray, optional
        Ice mixing ratio (kg/kg). Read from files if not given
    qs : 3-D ndarray, optional
        Snow mixing ratio (kg/kg). Read from files if not given
    qg : 3-D ndarray, optional
        Graupel mixing ratio (kg/kg). Read from files if not given
    t : int or <datetime.datetime> class or None, optional
        Time to read
        * None -- all times (defalut)
    dryrun : bool, optional
        * True -- dry run mode, only reading necessary data using 'sio' class
        * False -- do real computation (default)

    Returns
    -------
    qhydro : 3-D ndarray, optional
        Mixing ratio of all hydrometers (kg/kg)
    """
    qc_ = (sio.readvar('QC', t=t) if qc is None else qc)
    qr_ = (sio.readvar('QR', t=t) if qr is None else qr)
    qi_ = (sio.readvar('QI', t=t) if qi is None else qi)
    qs_ = (sio.readvar('QS', t=t) if qs is None else qs)
    qg_ = (sio.readvar('QG', t=t) if qg is None else qg)
    if dryrun:
        return None

    qhydro = qc_ + qr_ + qi_ + qs_ + qg_
    return qhydro


def calc_pt(sio, rho=None, rhot=None, qv=None, qhydro=None, tout=True, thetaout=False, t=None, dryrun=False):
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
    dryrun : bool, optional
        * True -- dry run mode, only reading necessary data using 'sio' class
        * False -- do real computation (default)

    Returns
    -------
    p : 3-D ndarray
        Pressure (Pa)
    t : 3-D ndarray, optional
        Temperature (K)
    theta : 3-D ndarray, optional
        Potential temperature (K)
    """
    rho_ = (sio.readvar('DENS', t=t) if rho is None else rho)
    rhot_ = (sio.readvar('RHOT', t=t) if rhot is None else rhot)
    qv_ = (sio.readvar('QV', t=t) if qv is None else qv)
    qhydro_ = (calc_qhydro(sio, t=t, dryrun=dryrun) if qhydro is None else qhydro)
    if dryrun:
        res = [None]
        if tout:
            res.append(None)
        if thetaout:
            res.append(None)
        return res

    Rdry = 287.04
    Rvap = 461.46
    CVdry = 717.60
    PRE00 = 100000.0

    Rtot = Rdry * (1. - qv_ - qhydro_) + Rvap * qv_
    CPovCV = ( CVdry + Rtot ) / CVdry

    p = PRE00 * np.power(rhot_ * Rtot / PRE00, CPovCV)
    res = [p]
    if tout:
        t = p / (rho_ * Rtot)
        res.append(t)
    if thetaout:
        theta = rhot_ / rho_
        res.append(theta)

    return res


def calc_ref(sio, min_dbz=-20., rho=None, qr=None, qs=None, qg=None, t=None, dryrun=False):
    """
    Calculate radar reflectivity

    Parameters
    ----------
    sio : <scale.io.ScaleIO> class
        Split SCALE I/O class
    min_dbz : float
        Minimum value of the reflectivity
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
    dryrun : bool, optional
        * True -- dry run mode, only reading necessary data using 'sio' class
        * False -- do real computation (default)

    Returns
    -------
    dbz : 3-D ndarray
        Radar reflectivity (dBZ)
    max_dbz : 2-D ndarray
        Maximum radar reflectivity in vertical (dBZ)
    """
    rho_ = (sio.readvar('DENS', t=t) if rho is None else rho)
    qr_ = (sio.readvar('QR', t=t) if qr is None else qr)
    qs_ = (sio.readvar('QS', t=t) if qs is None else qs)
    qg_ = (sio.readvar('QG', t=t) if qg is None else qg)
    if dryrun:
        return None, None

    qr_[qr_ < 1.e-10] = 1.e-10
    qs_[qs_ < 1.e-10] = 1.e-10
    qg_[qg_ < 1.e-10] = 1.e-10
    ref = 2.53e4 * (rho_ * qr_ * 1.0e3) ** 1.84 \
        + 3.48e3 * (rho_ * qs_ * 1.0e3) ** 1.66 \
        + 8.18e4 * (rho_ * qg_ * 1.0e3) ** 1.50
    dbz = 10. * np.log10(ref)
    dbz[dbz < min_dbz] = min_dbz
    max_dbz = ma.max(dbz, axis=0)

    return dbz, max_dbz


def extrap_z_t0(sio, temp, lprate=0.0065, zfree=200., height=None, t=None, dryrun=False):
    """
    Calculate smoothed lowest-level surface temperature extrapolated from 
    the free atmosphere

    Parameters
    ----------
    sio : <scale.io.ScaleIO> class
        Split SCALE I/O class
    temp : 3-D ndarray
        Input temperature (K)
    lprate : float, optional
        Assumed lapse rate (K/m). Default: 0.005
    zfree : float, optional
        Reference height of free atmosphere above the surface (m).
        Default: 1000.
    height : 3-D ndarray, optional
        Height in current grids (m). Read from files if not given
    t : int or <datetime.datetime> class or None, optional
        Time to read
        * None -- all times (defalut)
    dryrun : bool, optional
        * True -- dry run mode, only reading necessary data using 'sio' class
        * False -- do real computation (default)

    Returns
    -------
    t0_ext : 2-D ndarray
        Smoothed lowest-level temperature extrapolated from the free atmosphere (K)
    """
    height_ = (sio.readvar('height', t=t) if height is None else height)
    if dryrun:
        return None

#    height_[0,:,:] -= 1.e-2

    varshape = list(temp.shape)
    t0_ext = np.zeros(varshape[1:3], dtype=temp.dtype)

    for j in range(varshape[1]):
        for i in range(varshape[2]):
            t_ref = np.interp(height_[0,j,i]+zfree, np.copy(height_[:,j,i]), np.copy(temp[:,j,i]))
            t0_ext[j,i] = t_ref + lprate * zfree
    return t0_ext


def extrap_z_pt(sio, p0, t0_ext, height, qv=None, qhydro=None, lprate=0.0065, t=None, p=None, tk=None, theta=None, rho=None, rhot=None, dryrun=False):
    """
    Calculate extrapolated 3-D variables under the surface for z-level data

    Parameters
    ----------
    sio : <scale.io.ScaleIO> class
        Split SCALE I/O class
    p0 : 2-D ndarray
        ......
    t0_ext : 2-D ndarray
        ......
    height : 3-D ndarray
        ......
    qv : 3-D ndarray, optional
        Water vapor mixing ratio (kg/kg). Read from files if not given
    qhydro : 3-D ndarray, optional
        Mixing ratio of all hydrometers (kg/kg). Read from files if not given
    lprate : float, optional
        ......
    t : int or <datetime.datetime> class or None, optional
        Time to read
        * None -- all times (defalut)
    p : 3-D ndarray, optional, return
        Extrapolated pressure (Pa)
        * None -- do not calculate this variable
    tk : 3-D ndarray, optional, return
        Extrapolated temperature (K)
        * None -- do not calculate this variable
    theta : 3-D ndarray, optional, return
        Extrapolated potential temperature (K)
        * None -- do not calculate this variable
    rho : 3-D ndarray, optional, return
        Extrapolated density (kg/m3) (K)
        * None -- do not calculate this variable
    rhot : 3-D ndarray, optional, return
        Extrapolated rho * theta (kg/m3*K) (K)
        * None -- do not calculate this variable
    dryrun : bool, optional
        * True -- dry run mode, only reading necessary data using 'sio' class
        * False -- do real computation (default)
    """
    qv_ = (sio.readvar('QV', t=t) if qv is None else qv)
    qhydro_ = (calc_qhydro(sio, t=t, dryrun=dryrun) if qhydro is None else qhydro)
    if dryrun:
        return

    g = 9.80665
    Rdry = 287.04
    Rvap = 461.46
    CVdry = 717.60
    PRE00 = 100000.0

    Rtot = Rdry * (1. - qv_[0] - qhydro_[0]) + Rvap * qv_[0]
    RovCP = Rtot / (CVdry + Rtot)

    varshape = list(height.shape)

    for j in range(varshape[1]):
        for i in range(varshape[2]):
            for k in range(varshape[0]):
                if sio.z[k] <= height[0,j,i]:
                    tk_s = t0_ext[j,i] + lprate * (height[0,j,i] - sio.z[k])
                    p_s = p0[j,i] * (tk_s / t0_ext[j,i]) ** (g/Rtot[j,i]/lprate)
                    theta_s = tk_s * (PRE00 / p_s) ** RovCP[j,i]
                    rho_s = p_s / (Rtot[j,i] * tk_s)
                    if rhot is not None: rhot_s = rho_s * theta_s

                    if p     is not None: p[k,j,i] = p_s
                    if tk    is not None: tk[k,j,i] = tk_s
                    if theta is not None: theta[k,j,i] = theta_s
                    if rho   is not None: rho[k,j,i] = rho_s
                    if rhot  is not None: rhot[k,j,i] = rhot_s
#    return


def extrap_p_zt(sio, plevels, p, t0_ext, height, qv=None, qhydro=None, lprate=0.0065, t=None, z=None, tk=None, theta=None, rho=None, rhot=None, dryrun=False):
    """
    Calculate extrapolated 3-D variables under the surface for p-level data

    Parameters
    ----------
    sio : <scale.io.ScaleIO> class
        Split SCALE I/O class
    plevels : array_like
        ......
    p : 3-D ndarray
        ......
    t0_ext : 2-D ndarray
        ......
    height : 3-D ndarray
        ......
    qv : 3-D ndarray, optional
        Water vapor mixing ratio (kg/kg). Read from files if not given
    qhydro : 3-D ndarray, optional
        Mixing ratio of all hydrometers (kg/kg). Read from files if not given
    lprate : float, optional
        ......
    t : int or <datetime.datetime> class or None, optional
        Time to read
        * None -- all times (defalut)
    z : 3-D ndarray, optional, return
        Extrapolated height (m)
        * None -- do not calculate this variable
    tk : 3-D ndarray, optional, return
        Extrapolated temperature (K)
        * None -- do not calculate this variable
    theta : 3-D ndarray, optional, return
        Extrapolated potential temperature (K)
        * None -- do not calculate this variable
    rho : 3-D ndarray, optional, return
        Extrapolated density (kg/m3) (K)
        * None -- do not calculate this variable
    rhot : 3-D ndarray, optional, return
        Extrapolated rho * theta (kg/m3*K) (K)
        * None -- do not calculate this variable
    dryrun : bool, optional
        * True -- dry run mode, only reading necessary data using 'sio' class
        * False -- do real computation (default)
    """
    qv_ = (sio.readvar('QV', t=t) if qv is None else qv)
    qhydro_ = (calc_qhydro(sio, t=t, dryrun=dryrun) if qhydro is None else qhydro)
    if dryrun:
        return

    g = 9.80665
    Rdry = 287.04
    Rvap = 461.46
    CVdry = 717.60
    PRE00 = 100000.0

    Rtot = Rdry * (1. - qv_[0] - qhydro_[0]) + Rvap * qv_[0]
    RovCP = Rtot / (CVdry + Rtot)

    varshape = [len(plevels), height.shape[1], height.shape[2]]

    for j in range(varshape[1]):
        for i in range(varshape[2]):
            for k in range(varshape[0]):
                if plevels[k] >= p[0,j,i]:
                    if z is not None: z_s = height[0,j,i] + t0_ext[j,i]/lprate * (1. - (plevels[k]/p[0,j,i]) ** (Rtot[j,i]*lprate/g))
                    tk_s = t0_ext[j,i] * (plevels[k] / p[0,j,i]) ** (Rtot[j,i]*lprate/g)
                    theta_s = tk_s * (PRE00 / plevels[k]) ** RovCP[j,i]
                    rho_s = plevels[k] / (Rtot[j,i] * tk_s)
                    if rhot is not None: rhot_s = rho_s * theta_s

                    if z     is not None: z[k,j,i] = z_s
                    if tk    is not None: tk[k,j,i] = tk_s
                    if theta is not None: theta[k,j,i] = theta_s
                    if rho   is not None: rho[k,j,i] = rho_s
                    if rhot  is not None: rhot[k,j,i] = rhot_s
#    return


def calc_slp(sio, p0, t0_ext, height, qv=None, qhydro=None, lprate=0.0065, t=None, dryrun=False):
    """
    Calculate sea level pressure

    sio : <scale.io.ScaleIO> class
        Split SCALE I/O class
    p0 : 2-D ndarray
        ......
    t0_ext : 2-D ndarray
        ......
    height : 3-D ndarray
        ......
    qv : 3-D ndarray, optional
        Water vapor mixing ratio (kg/kg). Read from files if not given
    qhydro : 3-D ndarray, optional
        Mixing ratio of all hydrometers (kg/kg). Read from files if not given
    lprate : float, optional
        ......
    t : int or <datetime.datetime> class or None, optional
        Time to read
        * None -- all times (defalut)
    dryrun : bool, optional
        * True -- dry run mode, only reading necessary data using 'sio' class
        * False -- do real computation (default)

    Returns
    -------
    slp : 2-D ndarray
        Sea level pressure (Pa)
    """
    qv_ = (sio.readvar('QV', t=t) if qv is None else qv)
    qhydro_ = (calc_qhydro(sio, t=t, dryrun=dryrun) if qhydro is None else qhydro)
    if dryrun:
        return None

    g = 9.80665
    Rdry = 287.04
    Rvap = 461.46
    PRE00 = 100000.0

    Rtot = Rdry * (1. - qv_[0] - qhydro_[0]) + Rvap * qv_[0]

    varshape = list(qv_.shape)
    slp = np.zeros(varshape[1:3], dtype=qv_.dtype)

    for j in range(varshape[1]):
        for i in range(varshape[2]):
            slp[j,i] = p0[j,i] * (1 + lprate * height[0,j,i] / t0_ext[j,i]) ** (g/Rtot[j,i]/lprate)
    return slp


def lagrange_interp(x, xp, fp):
    """
    """
    res = ((x-xp[1]) * (x-xp[2])) / ((xp[0]-xp[1]) * (xp[0]-xp[2])) * fp[0] \
        + ((x-xp[0]) * (x-xp[2])) / ((xp[1]-xp[0]) * (xp[1]-xp[2])) * fp[1] \
        + ((x-xp[0]) * (x-xp[1])) / ((xp[2]-xp[0]) * (xp[2]-xp[1])) * fp[2]
    return res


def calc_rhosfc_psfc(sio, rho, pres, height, topo, t=None):
    """
    """
    g = 9.80665

    varshape = list(rho.shape)
    rhosfc = np.zeros(varshape[1:3], dtype=rho.dtype)
    psfc = np.zeros(varshape[1:3], dtype=rho.dtype)

    for j in range(varshape[1]):
        for i in range(varshape[2]):
            rhosfc[j,i] = lagrange_interp(topo[j,i], height[:,j,i], rho[:,j,i])
            psfc[j,i] = pres[0,j,i] + 0.5 * (rhosfc[j,i] + rho[0,j,i]) * g * (height[0,j,i] - topo[j,i])
    return rhosfc, psfc


#subroutine getgph (km, dp, tv, ps, orog, gph, pstag)
#!-------------------------------------------------------------------------------

#  implicit none

#  integer, intent(in)       :: km
#  real(r_sngl), intent(in)  :: dp(km), tv(km)
#  real(r_sngl), intent(in)  :: ps, orog
#  real(r_sngl), intent(out) :: gph(km+1)
#  real(r_sngl), intent(out) :: pstag(km+1)
#  integer :: k

#!-------------------------------------------------------------------------------

#  gph(1) = orog
#  pstag(1) = ps
#  do k = 1, km-1
#    pstag(k+1) = pstag(k) - dp(k)
#    gph(k+1) = gph(k) + con_rog * tv(k) * log(pstag(k) / pstag(k+1))
#  enddo
#  pstag(km+1) = 0.
#  gph(km+1) = 9.99e33

#!-------------------------------------------------------------------------------
#end subroutine getgph



#subroutine getrh (km, p, sh, t, rh)
#!-------------------------------------------------------------------------------

#  implicit none

#  integer, intent(in)  :: km
#  real(r_sngl), intent(in)  :: p(km), sh(km), t(km)
#  real(r_sngl), intent(out) :: rh(km)
#  real(r_sngl) :: tr, es, shs
#  integer :: k

#!-------------------------------------------------------------------------------

#  do k = 1, km
#    tr = con_ttp / t(k)
#    es = con_psat * (tr**con_xpona) * exp(con_xponb*(1.-tr))
#    es = min(es, p(k))
#    shs = con_eps * es / (p(k) + con_epsm1 * es)
#    rh(k) = 1.e2 * min(max(sh(k)/shs, 0.), 1.)
#!    rh(k) = 1.e2 * sh(k) / shs
#  enddo

#!-------------------------------------------------------------------------------
#end subroutine getrh
