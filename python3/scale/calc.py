import numpy as np
import numpy.ma as ma
import numpy.testing as npt
from mpl_toolkits.basemap import Basemap


__all__ = ['set_bmap', 'calc_rotate_winds', 
           'calc_height', 'interp_z', 'interp_p', 'calc_destagger', 'calc_destagger_uvw', 'calc_qhydro', 'calc_pt',
           'calc_ref', 'extrap_z_t0', 'extrap_z_pt', 'extrap_p_zt', 'calc_slp', 'calc_rhosfc_psfc']


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


def calc_rotate_winds(sio, bmap, u=None, v=None, t=None):
    """
    Calculate the rotation of u, v winds

    XXXXXX
    """
    if u is None:
        u = sio.readvar('U', t=t)
    if v is None:
        v = sio.readvar('V', t=t)

    ny, nx = sio.lon.shape
    nz = len(sio.z)
    lon = sio.lon[sio.bufsize:ny-sio.bufsize, sio.bufsize:nx-sio.bufsize]
    lat = sio.lon[sio.bufsize:ny-sio.bufsize, sio.bufsize:nx-sio.bufsize]
    if nz == 3:
        lon = np.repeat(lon[np.newaxis,:,:], nz, axis=0)
        lat = np.repeat(lon[np.newaxis,:,:], nz, axis=0)

    u_rot, v_rot = bmap.rotate_vector(u3d.ravel(), v3d.ravel(), lon.ravel(), lat.ravel())
    if nz == 3:
        u_rot = u_rot.reshape(nz, ny, nx)
        v_rot = v_rot.reshape(nz, ny, nx)
    else:
        u_rot = u_rot.reshape(ny, nx)
        v_rot = v_rot.reshape(ny, nx)

    mag_square = u * u + v * v
    tmpcos = (u_rot * u + v_rot * v) / mag_square
    tmpsin = (u_rot * v - v_rot * u) / mag_square
    u_rot = tmpcos * u - tmpsin * v
    v_rot = tmpsin * u + tmpcos * v

    return u_rot, v_rot


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
        topo = sio.readvar('TOPO')

    height = np.zeros((sio.dimdef['len']['z'][0], topo.shape[0], topo.shape[1]), dtype=topo.dtype)
    height_h = np.zeros((sio.dimdef['len']['zh'][0], topo.shape[0], topo.shape[1]), dtype=topo.dtype)

    for k in range(len(sio.z)):
        height[k,:,:] = topo + (sio.zh[-1] - topo) / sio.zh[-1] * sio.z[k]
    for k in range(len(sio.zh)):
        height_h[k,:,:] = topo + (sio.zh[-1] - topo) / sio.zh[-1] * sio.zh[k]

    return height, height_h


def interp_z(sio, var, height=None, t=None, extrap=False):
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

    Returns
    -------
    varout : 3-D ndarray
        Interpolated variable
    """
    if height is None:
        height = sio.readvar('height', t=t)
    height[0,:,:] -= 1.e-2

    varshape = list(var.shape)
    varshape[0] = sio.dimdef['len']['z'][0]
    varout = ma.masked_all(varshape, dtype=var.dtype)
    if type(var) == ma.MaskedArray:
        varout.fill_value = var.fill_value

    for j in range(varshape[1]):
        for i in range(varshape[2]):
            if extrap:
                varout[:,j,i] = np.interp(sio.z, np.copy(height[:,j,i]), np.copy(var[:,j,i]), right=-9.99e+33)
#                varout[:,j,i] = np.interp(sio.z, height[:,j,i], var[:,j,i], right=-9.99e+33)
            else:
                varout[:,j,i] = np.interp(sio.z, np.copy(height[:,j,i]), np.copy(var[:,j,i]), left=-9.99e+33, right=-9.99e+33)
    varout.mask[varout == -9.99e+33] = True

    return varout


#def interp_p_thread(sequence):
#    if sequence[3]:
#        return np.interp(sequence[0], sequence[1], sequence[2], left=-9.99e+33)
#    else:
#        return np.interp(sequence[0], sequence[1], sequence[2], left=-9.99e+33, right=-9.99e+33)


def interp_p(sio, var, plevels, p=None, t=None, extrap=False):
#def interp_p(sio, var, plevels, p=None, t=None, extrap=False, threads=1):
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
    extrap : bool, optional
        Extrapolate low values?
        * True -- do extrapolation
        * False -- do not do extrapolation (default)

    Returns
    -------
    varout : 3-D ndarray
        Interpolated variable
    """

    if p is None:
        p = calc_pt(sio, tout=False, thetaout=False, t=t)[0]
    varshape = list(var.shape)
    varshape[0] = len(plevels)
    varout = ma.masked_all(varshape, dtype=var.dtype)
    if type(var) == ma.MaskedArray:
        varout.fill_value = var.fill_value

    log_plevels = np.log(plevels)
    log_p_inv = np.log(np.transpose(p, (1, 2, 0))[:,:,::-1])
    var_inv = np.copy(np.transpose(var, (1, 2, 0))[:,:,::-1])

#    if threads == 1:
    for j in range(varshape[1]):
        for i in range(varshape[2]):
            if extrap:
                varout[:,j,i] = np.interp(log_plevels, log_p_inv[j,i,:], var_inv[j,i,:], left=-9.99e+33)
            else:
                varout[:,j,i] = np.interp(log_plevels, log_p_inv[j,i,:], var_inv[j,i,:], left=-9.99e+33, right=-9.99e+33)
#    else:
##        def tmpfunc(ij):
###            i = ij % varshape[2]
###            j = ij // varshape[2]
##            return [1., 2., 4.]
##            if extrap:
##                return np.interp(log_plevels, log_p_inv[j,i,:], var_inv[j,i,:], left=-9.99e+33)
##            else:
##                return np.interp(log_plevels, log_p_inv[j,i,:], var_inv[j,i,:], left=-9.99e+33, right=-9.99e+33)

#        sequences = []
#        for j in range(varshape[1]):
#            for i in range(varshape[2]):
#                sequences.append((log_plevels, log_p_inv[j,i,:], var_inv[j,i,:], extrap))

#        pool = Pool(processes=threads)
#        results = pool.map(interp_p_thread, sequences) #, chunksize=100)
##        cleaned = [x for x in results if not x is None]
## not optimal but safe
#        pool.close()
#        pool.join()

#        for j in range(varshape[1]):
#            for i in range(varshape[2]):
#                ij = j * varshape[2] + i
#                varout[:,j,i] = results[ij]


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


def calc_destagger_uvw(sio, rho=None, momx=None, momy=None, momz=None, destagger=True, first_grd=True, t=None):
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
    if rho is None:
        rho = sio.readvar('DENS', t=t)
    if momx is None:
        if first_grd:
            momx = sio.readvar('MOMX', t=t, bufsize=0)[:,sio.bufsize:-sio.bufsize,sio.bufsize-1:-sio.bufsize]
        else:
            momx = sio.readvar('MOMX', t=t)
    if momy is None:
        if first_grd:
            momy = sio.readvar('MOMY', t=t, bufsize=0)[:,sio.bufsize-1:-sio.bufsize,sio.bufsize:-sio.bufsize]
        else:
            momy = sio.readvar('MOMY', t=t)
    if momz is None:
        momz = sio.readvar('MOMZ', t=t)

    if destagger:
#        print(' --- destagger momx')
        momx = calc_destagger(momx, axis=2, first_grd=first_grd)
#        print(' --- destagger momy')
        momy = calc_destagger(momy, axis=1, first_grd=first_grd)
#        print(' --- destagger momz')
        momz = calc_destagger(momz, axis=0, first_grd=False)
    momz[0,:,:] = 0.

    u = momx / rho
    v = momy / rho
    w = momz / rho

    return u, v, w, momx, momy, momz


def calc_qhydro(sio, qc=None, qr=None, qi=None, qs=None, qg=None, t=None):
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

    Returns
    -------
    qhydro : 3-D ndarray, optional
        Mixing ratio of all hydrometers (kg/kg)
    """
    if qc is None:
        qc = sio.readvar('QC', t=t)
    if qr is None:
        qr = sio.readvar('QR', t=t)
    if qi is None:
        qi = sio.readvar('QI', t=t)
    if qs is None:
        qs = sio.readvar('QS', t=t)
    if qg is None:
        qg = sio.readvar('QG', t=t)

    qhydro = qc + qr + qi + qs + qg
    return qhydro


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
        rho = sio.readvar('DENS', t=t)
    if rhot is None:
        rhot = sio.readvar('RHOT', t=t)
    if qv is None:
        qv = sio.readvar('QV', t=t)
    if qhydro is None:
        qhydro = calc_qhydro(sio, t=t)

    Rdry = 287.04
    Rvap = 461.46
    CVdry = 717.60
    PRE00 = 100000.0

    Rtot = Rdry * (1. - qv - qhydro) + Rvap * qv
    CPovCV = ( CVdry + Rtot ) / CVdry

    p = PRE00 * np.power(rhot * Rtot / PRE00, CPovCV)
    res = [p]
    if tout:
        t = p / (rho * Rtot)
        res.append(t)
    if thetaout:
        theta = rhot / rho
        res.append(theta)

    return res


def calc_ref(sio, min_dbz=-20., rho=None, qr=None, qs=None, qg=None, t=None):
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

    Returns
    -------
    dbz : 3-D ndarray
        Radar reflectivity (dBZ)
    max_dbz : 2-D ndarray
        Maximum radar reflectivity in vertical (dBZ)
    """
    if rho is None:
        rho = sio.readvar('DENS', t=t)
    if qr is None:
        qr = sio.readvar('QR', t=t)
    if qs is None:
        qs = sio.readvar('QS', t=t)
    if qg is None:
        qg = sio.readvar('QG', t=t)

    qr[qr < 1.e-10] = 1.e-10
    qs[qs < 1.e-10] = 1.e-10
    qg[qg < 1.e-10] = 1.e-10
    ref = 2.53e4 * (rho * qr * 1.0e3) ** 1.84 \
        + 3.48e3 * (rho * qs * 1.0e3) ** 1.66 \
        + 8.18e4 * (rho * qg * 1.0e3) ** 1.50
    dbz = 10. * np.log10(ref)
    dbz[dbz < min_dbz] = min_dbz
    max_dbz = ma.max(dbz, axis=0)

    return dbz, max_dbz


def extrap_z_t0(sio, temp, lprate=0.005, zfree=1000., height=None, t=None):
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
    zfree : float,optional
        Reference height of free atmosphere above the surface (m).
        Default: 1000.
    height : 3-D ndarray, optional
        Height in current grids (m). Read from files if not given
    t : int or <datetime.datetime> class or None, optional
        Time to read
        * None -- all times (defalut)

    Returns
    -------
    t0_ext : 2-D ndarray
        Smoothed lowest-level temperature extrapolated from the free atmosphere (K)
    """
    if height is None:
        height = sio.readvar('height', t=t)
#    height[0,:,:] -= 1.e-2

    varshape = list(temp.shape)
    t0_ext = np.zeros(varshape[1:3], dtype=temp.dtype)

    for j in range(varshape[1]):
        for i in range(varshape[2]):
            t_ref = np.interp(height[0,j,i]+zfree, np.copy(height[:,j,i]), np.copy(temp[:,j,i]))
            t0_ext[j,i] = t_ref + lprate * zfree
    return t0_ext


def extrap_z_pt(sio, qv, qhydro, p0, t0_ext, height, lprate=0.005, t=None, p=None, tk=None, theta=None, rho=None, rhot=None):
#def extrap_z_pt(sio, qv=None, qhydro=None, p0, t0_ext, height, lprate=0.005, t=None, p=None, tk=None, theta=None, rho=None, rhot=None):
    """
    Calculate extrapolated 3-D variables under the surface for z-level data

    Parameters
    ----------
    sio : <scale.io.ScaleIO> class
        Split SCALE I/O class
    qv : 3-D ndarray, optional
        Water vapor mixing ratio (kg/kg). Read from files if not given
    qhydro : 3-D ndarray, optional
        Mixing ratio of all hydrometers (kg/kg). Read from files if not given
    p0 : 2-D ndarray
        ......
    t0_ext : 2-D ndarray
        ......
    height : 3-D ndarray
        ......
    lprate : float
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
    """
    if qv is None:
        qv = sio.readvar('QV', t=t)
    if qhydro is None:
        qhydro = calc_qhydro(sio, t=t)

    g = 9.80665
    Rdry = 287.04
    Rvap = 461.46
    CVdry = 717.60
    PRE00 = 100000.0

    Rtot = Rdry * (1. - qv[0] - qhydro[0]) + Rvap * qv[0]
    RovCP = Rtot / (CVdry + Rtot)

    varshape = list(p.shape)

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


def extrap_p_zt(sio, plevels, qv, qhydro, p, t0_ext, height, lprate=0.005, t=None, z=None, tk=None, theta=None, rho=None, rhot=None):
#def extrap_p_zt(sio, plevels, qv=None, qhydro=None, p, t0_ext, height, lprate=0.005, t=None, z=None, tk=None, theta=None, rho=None, rhot=None):
    """
    Calculate extrapolated 3-D variables under the surface for p-level data

    Parameters
    ----------
    sio : <scale.io.ScaleIO> class
        Split SCALE I/O class
    plevels : array_like
        ......
    qv : 3-D ndarray, optional
        Water vapor mixing ratio (kg/kg). Read from files if not given
    qhydro : 3-D ndarray, optional
        Mixing ratio of all hydrometers (kg/kg). Read from files if not given
    p : 3-D ndarray
        ......
    t0_ext : 2-D ndarray
        ......
    height : 3-D ndarray
        ......
    lprate : float
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
    """
    if qv is None:
        qv = sio.readvar('QV', t=t)
    if qhydro is None:
        qhydro = calc_qhydro(sio, t=t)

    g = 9.80665
    Rdry = 287.04
    Rvap = 461.46
    CVdry = 717.60
    PRE00 = 100000.0

    Rtot = Rdry * (1. - qv[0] - qhydro[0]) + Rvap * qv[0]
    RovCP = Rtot / (CVdry + Rtot)

    varshape = list(z.shape)

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


def calc_slp(sio, qv, qhydro, p0, t0_ext, height, lprate=0.005, t=None):
#def calc_slp(sio, qv=None, qhydro=None, p0, t0_ext, height, lprate=0.005, t=None):
    """
    Calculate sea level pressure

    sio : <scale.io.ScaleIO> class
        Split SCALE I/O class
    qv : 3-D ndarray, optional
        Water vapor mixing ratio (kg/kg). Read from files if not given
    qhydro : 3-D ndarray, optional
        Mixing ratio of all hydrometers (kg/kg). Read from files if not given
    p0 : 2-D ndarray
        ......
    t0_ext : 2-D ndarray
        ......
    height : 3-D ndarray
        ......
    lprate : float
        ......
    t : int or <datetime.datetime> class or None, optional
        Time to read
        * None -- all times (defalut)

    Returns
    -------
    slp : 2-D ndarray
        Sea level pressure (Pa)
    """
    if qv is None:
        qv = sio.readvar('QV', t=t)
    if qhydro is None:
        qhydro = calc_qhydro(sio, t=t)

    g = 9.80665
    Rdry = 287.04
    Rvap = 461.46
    PRE00 = 100000.0

    Rtot = Rdry * (1. - qv[0] - qhydro[0]) + Rvap * qv[0]

    varshape = list(qv.shape)
    slp = np.zeros(varshape[1:3], dtype=qv.dtype)

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
