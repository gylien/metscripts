import numpy as np
import numpy.ma as ma
import datetime as dt
import os
import time
from mpl_toolkits.basemap import Basemap
from scipy.interpolate import interpn
from .io import ScaleIO
from .calc import *
#from scale.proj import *
import gradsio
import sys

#try:
#    from mpi4py import MPI
#except ImportError:
#    pass


__all__ = ['convert_hintp', 'convert']


conf_default = {
'ftype': 'restart',
'missing': -9.99e33,
'vcoor': 'p',
'plevels': [100000., 92500., 85000., 70000., 50000., 30000., 20000., 10000., 5000., 2000., 1000.],
'varout_3d': ['u', 'v', 'w', 'p', 'tk', 'theta', 'rho', 'momx', 'momy', 'momz', 'rhot', 'z', 'qv', 'qc', 'qr', 'qi', 'qs', 'qg', 'qhydro', 'dbz'],
'varout_2d': ['topo', 'rhosfc', 'psfc', 'slp', 'rain', 'snow', 'max_dbz', 'glon', 'glat'],
'proj': {'type': 'LC',
         'basepoint_lon': 135.,
         'basepoint_lat': 35.,
         'basepoint_x': None,
         'basepoint_y': None,
         'LC_lat1': 30.,
         'LC_lat2': 40.,},
'extrap': True,
'lprate': 0.0065,
'zfree': 200.,
'windrot': True,
'dlon': 1.,
'dlat': 1.,
'tstart': 0,
'tend': -1,
'tskip': 1,
}

var_3d_name = {
'u': 'u-wind (m/s)',
'v': 'v-wind (m/s)',
'w': 'w-wind (m/s)',
'p': 'Pressure (Pa)',
'tk': 'Temperature (K)',
'theta': 'Potential temperature (K)',
'rho': 'Air density (kg/m^3)',
'momx': 'x-momentum (kg/m2/s)',
'momy': 'y-momentum (kg/m2/s)',
'momz': 'z-momentum (kg/m2/s)',
'rhot': 'rho * theta (kg/m3*K)',
'z': 'Height (m)',
'qv': 'Water vapor mixing ratio (kg/kg)',
'qc': 'Cloud water mixing ratio (kg/kg)',
'qr': 'Rain water mixing ratio (kg/kg)',
'qi': 'Cloud ice mixing ratio (kg/kg)',
'qs': 'Snow mixing ratio (kg/kg)',
'qg': 'Graupel mixing ratio (kg/kg)',
'qhydro': 'Mixing ratio of all hydrometers (kg/kg)',
'rh': 'Relative humidity (%)',
'dbz': 'Radar reflectivity (dBZ)'
}
var_2d_name = {
'topo': 'Topography height (m)',
'u10': '10m u-wind (m/s)',
'v10': '10m v-wind (m/s)',
't2': '2m temperature (K)',
'q2': '2m water vapor mixing ratio (kg/kg)',
'rhosfc': 'Surface air density (kg/m^3)',
'psfc': 'Surface pressure (Pa)',
'slp': 'Sea level pressure (Pa)',
'rain': 'Surface rain rate (mm/s)',
'snow': 'Surface snow rate (mm/s)',
'max_dbz': 'Maximum radar reflectivity (dBZ)',
'olr': 'TOA net longwave radiation flux (W/m^2)',
'tsfc': 'Surface skin temperature (merged) (K)',
'tsfcocean': 'Ocean surface skin temperature (K)',
'sst': 'Temperature at uppermost ocean layer (K)',
'glon': 'Longitude (degree)',
'glat': 'Latitude (degree)'
}

var_3d = list(var_3d_name.keys())
var_2d = list(var_2d_name.keys())
var_3d_sprd = ['u', 'v', 'w', 'tk', 'p', 'qv', 'qc', 'qr', 'qi', 'qs', 'qg']
var_2d_sprd = []


def rc(**kwargs):
    """
    """
    from copy import deepcopy
    conf = deepcopy(conf_default)
    for key, value in list(kwargs.items()):
        if key in conf:
            if value == None:
                pass
            elif key == 'proj':
                for key2, value2 in list(value.items()):
                    if key2 in conf[key]:
                        conf[key][key2] = value2
                    else:
                        raise KeyError("'{0:s}' is not a key of conf['{1:s}'].".format(key2, key))
            else:
                conf[key] = value
        else:
            raise KeyError("'{0:s}' is not a configuration key.".format(key))
    return conf


def convert_readvar(sio, bmap, topo, conf, var_necessary, it, tskip_a, dryrun=False, myrankmsg=''):
    """
    """
    X = {}
    t0_ext = None
    X['topo'] = topo

    if conf['ftype'] == 'restart':
        for ivar, ivarf in ('rho', 'DENS'), ('rhot', 'RHOT'), \
                           ('qv', 'QV'), ('qc', 'QC'), ('qr', 'QR'), ('qi', 'QI'), ('qs', 'QS'), ('qg', 'QG'):
            if ivar in conf['varout_3d'] or ivar in var_necessary:
                X[ivar] = sio.readvar(ivarf, t=it)
        for ivar, ivarf in ('rain', 'SFLX_rain'), ('snow', 'SFLX_snow'), ('glon', 'lon'), ('glat', 'lat'):
            if ivar in conf['varout_2d'] or ivar in var_necessary:
                X[ivar] = sio.readvar(ivarf, t=it)

        if 'u' in conf['varout_3d'] or 'v' in conf['varout_3d']:
            if not dryrun:
                print(myrankmsg, 'Calculate: destaggered u, v, w, momx, momy, momz')
            X['u'], X['v'], X['w'], X['momx'], X['momy'], X['momz'] = \
                calc_destagger_uvw(sio, first_grd=True, t=it, dryrun=dryrun)
            if not dryrun:
                print(myrankmsg, 'Calculate: rotate u, v')
                X['u'], X['v'] = calc_rotate_winds(sio, bmap, u=X['u'], v=X['v'], t=it)

        if not dryrun:
            print(myrankmsg, 'Calculate: qhydro')
        X['qhydro'] = calc_qhydro(sio, t=it, dryrun=dryrun)

        if not dryrun:
            print(myrankmsg, 'Calculate: p, t, theta')
        X['p'], X['tk'], X['theta'] = calc_pt(sio, qhydro=X['qhydro'], tout=True, thetaout=True, t=it, dryrun=dryrun)

        if 'dbz' in conf['varout_3d'] or 'max_dbz' in conf['varout_2d']:
            if not dryrun:
                print(myrankmsg, 'Calculate: dbz, max_dbz')
            X['dbz'], X['max_dbz'] = calc_ref(sio, t=it, dryrun=dryrun)

        if not dryrun:
            print(myrankmsg, 'Calculate: z')
        X['z'], height_h = calc_height(sio, topo=X['topo'], dryrun=dryrun)

        if 'rhosfc' in conf['varout_2d'] or 'psfc' in conf['varout_2d']:
            if not dryrun:
                print(myrankmsg, 'Calculate: rhosfc, psfc')
                X['rhosfc'], X['psfc'] = calc_rhosfc_psfc(sio, rho=X['rho'], pres=X['p'], height=X['z'], topo=X['topo'], t=it)

        if 'slp' in conf['varout_2d'] or (conf['extrap'] and (conf['vcoor'] == 'z' or conf['vcoor'] == 'p')):
            if not dryrun:
                print(myrankmsg, 'Calculate smoothed lowest-level surface temperature extrapolated from the free atmosphere')
                t0_ext = extrap_z_t0(sio, X['tk'], lprate=conf['lprate'], zfree=conf['zfree'], height=X['z'], t=it)

            if 'slp' in conf['varout_2d']:
                if not dryrun:
                    print(myrankmsg, 'Calculate: slp')
                    X['slp'] = calc_slp(sio, p0=X['p'][0], t0_ext=t0_ext, height=X['z'], 
                                        qv=X['qv'], qhydro=X['qhydro'], lprate=conf['lprate'], t=it)

    elif conf['ftype'] == 'restart_sprd':
        for ivar, ivarf in ('u', 'DENS'), ('v', 'MOMX'), ('w', 'MOMY'), ('tk', 'MOMZ'), ('p', 'RHOT'), \
                           ('qv', 'QV'), ('qc', 'QC'), ('qr', 'QR'), ('qi', 'QI'), ('qs', 'QS'), ('qg', 'QG'):
            if ivar in conf['varout_3d']:
                X[ivar] = sio.readvar(ivarf, t=it)

        if not dryrun:
            print(myrankmsg, 'Destagger: u, v, w')
            if 'u' in conf['varout_3d']:
                X['u'] = calc_destagger(X['u'], axis=2, first_grd=False)
            if 'v' in conf['varout_3d']:
                X['v'] = calc_destagger(X['v'], axis=1, first_grd=False)
            if 'w' in conf['varout_3d']:
                X['w'] = calc_destagger(X['w'], axis=0, first_grd=False)
                X['w'][0,:,:] = 0.

    elif conf['ftype'] == 'history':
        run_calc_pt = False
        for ivar, ivarf in ('rho', 'DENS'), ('momx', 'MOMX'), ('momy', 'MOMY'), ('momz', 'MOMZ'), ('rhot', 'RHOT'), \
                           ('qv', 'QV'), ('qc', 'QC'), ('qr', 'QR'), ('qi', 'QI'), ('qs', 'QS'), ('qg', 'QG'), ('qhydro', 'QHYD'), \
                           ('u', 'U'), ('v', 'V'), ('w', 'W'), ('tk', 'T'), ('p', 'PRES'), ('theta', 'PT'), ('rh', 'RH'):
            if ivar in conf['varout_3d'] or ivar in var_necessary:
                try:
                    X[ivar] = sio.readvar(ivarf, t=it)
                except IOError as err:
                    if ivar == 'theta':
                        run_calc_pt = True
                        pass
                    else:
                        raise err

        if 'u' in conf['varout_3d'] or 'v' in conf['varout_3d']:
            if not dryrun:
                print(myrankmsg, 'Calculate: rotate u, v')
                X['u'], X['v'] = calc_rotate_winds(sio, bmap, u=X['u'], v=X['v'], t=it)

        if run_calc_pt:
            if not dryrun:
                print(myrankmsg, 'Calculate: p, t, theta')
            X['p'], X['tk'], X['theta'] = calc_pt(sio, qhydro=X['qhydro'], tout=True, thetaout=True, t=it, dryrun=dryrun)

        for ivar, ivarf in ('u10', 'U10'), ('v10', 'V10'), ('t2', 'T2'), ('q2', 'Q2'), ('olr', 'OLR'), ('slp', 'MSLP'), \
                           ('sst', 'OCEAN_TEMP'), ('tsfc', 'SFC_TEMP'), ('tsfcocean', 'OCEAN_SFC_TEMP'), ('glon', 'lon'), ('glat', 'lat'):
            if ivar in conf['varout_2d'] or ivar in var_necessary:
                X[ivar] = sio.readvar(ivarf, t=it)

        if 'u10' in conf['varout_2d'] or 'v10' in conf['varout_2d']:
            if not dryrun:
                print(myrankmsg, 'Calculate: rotate u10, v10')
                X['u10'], X['v10'] = calc_rotate_winds(sio, bmap, u=X['u10'], v=X['v10'], t=it)

        for ivar, ivarf in ('rain', 'RAIN'), ('snow', 'SNOW'):
            if ivar in conf['varout_2d'] or ivar in var_necessary:
                iits = max(it-tskip_a+1, 0)
                if dryrun:
                    for iit in range(iits, it+1):
                        sio.readvar(ivarf, t=iit)
                else:
                    X[ivar] = sio.readvar(ivarf, t=iits)
                    for iit in range(iits+1, it+1):
                        X[ivar] += sio.readvar(ivarf, t=iit)
                    X[ivar] /= (it - iits + 1)

        if 'qhydro' in conf['varout_3d'] and 'qhydro' not in X:
            if not dryrun:
                print(myrankmsg, 'Calculate: qhydro')
            X['qhydro'] = calc_qhydro(sio, t=it, dryrun=dryrun)
        if 'dbz' in conf['varout_3d'] or 'max_dbz' in conf['varout_2d']:
            if not dryrun:
                print(myrankmsg, 'Calculate: dbz, max_dbz')
            X['dbz'], X['max_dbz'] = calc_ref(sio, t=it, dryrun=dryrun)

        if not dryrun:
            print(myrankmsg, 'Calculate: z')
        X['z'], height_h = calc_height(sio, topo=X['topo'], dryrun=dryrun)
        if not dryrun:
            X['z'] = X['z'].astype('f4')

        if 'slp' in conf['varout_2d'] and 'slp' not in X or (conf['extrap'] and (conf['vcoor'] == 'z' or conf['vcoor'] == 'p')):
            if not dryrun:
                print(myrankmsg, 'Calculate smoothed lowest-level surface temperature extrapolated from the free atmosphere')
            t0_ext = extrap_z_t0(sio, X['tk'], lprate=conf['lprate'], zfree=conf['zfree'], height=X['z'], t=it, dryrun=dryrun)

            if 'slp' in conf['varout_2d'] and 'slp' not in X:
                if not dryrun:
                    print(myrankmsg, 'Calculate: slp')
                X['slp'] = calc_slp(sio, p0=X['p'][0], t0_ext=t0_ext, height=X['z'],
                                    qv=X['qv'], qhydro=X['qhydro'], lprate=conf['lprate'], t=it, dryrun=dryrun)

    elif conf['ftype'] == 'history_z':
        sys.exit('not done yet...')

    else:
        raise ValueError("ftype = '{0:s}' is not supported. ftype: {'restart', 'restart_sprd', 'history', 'history_z'}".format(conf['ftype']))

    return X, t0_ext


def convert_vintp(sio, conf, X, t0_ext, nx, ny, nzout, it, myrankmsg=''):
    """
    """
    Xitp = {}

    # do not consider 'history_z'...
    # require: X['z']
    #  <extrap> X['p'], X['tk'], X['theta'], X['rho'], X['rhot'], X['qv'], X['qhydro']
    if conf['vcoor'] == 'z' and conf['ftype'] != 'restart_sprd':
        for ivar in X:
#                    if ivar != 'z' and (ivar in conf['varout_3d'] or ivar in ['p', 'tk', 'theta', 'rho', 'rhot']):
            if ivar != 'z' and (ivar in conf['varout_3d']):
                print(myrankmsg, 'Vertical interpolation at Z-coordinate: ', ivar)
                Xitp[ivar] = interp_z(sio, X[ivar], height=X['z'], t=it, extrap=conf['extrap'])
        if 'z' in conf['varout_3d']:
            Xitp['z'] = np.empty((nzout, ny, nx), dtype=sio.z.dtype)
            for ilev in range(nzout):
                Xitp['z'][ilev] = sio.z[ilev]

        if conf['extrap']:
            kws = {}
            kwslist = ''
            for ivar in ['p', 'tk', 'theta', 'rho', 'rhot']:
                if ivar in conf['varout_3d']:
                    kws[ivar] = Xitp[ivar]
                    kwslist += ivar + ', '
            print(myrankmsg, ' Calculate extrapolated values under the surface assuming a constant lapse rate: ' + kwslist[0:-2])
            extrap_z_pt(sio, p0=X['p'][0], t0_ext=t0_ext, height=X['z'],
                        qv=X['qv'], qhydro=X['qhydro'], lprate=conf['lprate'], t=it, **kws)


    # do not consider 'history_z'...
    # require: X['p']
    #  <extrap> X['z'], X['tk'], X['theta'], X['rho'], X['rhot'], X['qv'], X['qhydro']
    if conf['vcoor'] == 'p' and conf['ftype'] != 'restart_sprd':

        for ivar in X:
            if ivar != 'p' and (ivar in conf['varout_3d']):
                print(myrankmsg, 'Vertical interpolation at P-coordinate: ', ivar)
                Xitp[ivar] = interp_p(sio, X[ivar], conf['plevels'], p=X['p'], t=it, extrap=conf['extrap'])

        if 'p' in conf['varout_3d']:
            varshape = list(X[conf['varout_3d'][0]].shape)
            varshape[0] = len(conf['plevels'])
            Xitp['p'] = np.empty(varshape, dtype=X[ivar].dtype)
            for ilev in range(len(conf['plevels'])):
                Xitp['p'][ilev] = conf['plevels'][ilev]

        if conf['extrap']:
            kws = {}
            kwslist = ''
            for ivar in ['z', 'tk', 'theta', 'rho', 'rhot']:
                if ivar in conf['varout_3d']:
                    kws[ivar] = Xitp[ivar]
                    kwslist += ivar + ', '
            print(myrankmsg, 'Calculate extrapolated values under the surface assuming a constant lapse rate: ' + kwslist[0:-2])
            extrap_p_zt(sio, conf['plevels'], p=X['p'], t0_ext=t0_ext, height=X['z'],
                        qv=X['qv'], qhydro=X['qhydro'], lprate=conf['lprate'], t=it, **kws)

    # Convert dictionary of numpy nd array to two grand numpy ndarrays for 3D and 2D data
    if conf['vcoor'] == 'o' or conf['ftype'] == 'restart_sprd':
        Xout = X
    else:
        Xout = Xitp
    nv3d = 0
    for ivar in conf['varout_3d']:
        if ivar not in Xout:
            raise ValueError("Output variable '" + ivar + "' has not been calculated.")
        nv3d += 1
    X3d = np.empty((nv3d, nzout, ny, nx), dtype='f4')
    iv3d = 0
    for ivar in conf['varout_3d']:
        if type(Xout[ivar]) == ma.MaskedArray:
            X3d[iv3d] = Xout[ivar].filled(fill_value=conf['missing'])
        else:
            X3d[iv3d] = Xout[ivar]
        iv3d += 1

    # always look for 'X' (instead of 'Xout') for 2D variables
    nv2d = 0
    for ivar in conf['varout_2d']:
        if ivar not in X:
            raise ValueError("Output variable '" + ivar + "' has not been calculated.")
        nv2d += 1
    X2d = np.empty((nv2d, ny, nx), dtype='f4')
    iv2d = 0
    for ivar in conf['varout_2d']:
        if type(X[ivar]) == ma.MaskedArray:
            X2d[iv2d] = X[ivar].filled(fill_value=conf['missing'])
        else:
            X2d[iv2d] = X[ivar]
        iv2d += 1

    return X3d, X2d, nv3d, nv2d


def convert_hintp(sio, bmap, conf, X3d, X2d, dlon, dlat, missing):
    """
    """
    nv3d = X3d.shape[0]
    nv2d = X2d.shape[0]
    nz = X3d.shape[1]
    ny = X3d.shape[2]
    nx = X3d.shape[3]

    dx = sio.dimdef['coor_g']['x'][sio.bufsize+1] - sio.dimdef['coor_g']['x'][sio.bufsize]
    dy = sio.dimdef['coor_g']['y'][sio.bufsize+1] - sio.dimdef['coor_g']['y'][sio.bufsize]

    lon = sio.readvar('lon')
    lat = sio.readvar('lat')
    lon_s = np.floor(np.min(lon) / dlon) * dlon
    lon_e = np.ceil(np.max(lon) / dlon) * dlon
    lat_s = np.floor(np.min(lat) / dlat) * dlat
    lat_e = np.ceil(np.max(lat) / dlat) * dlat

    lono1d = np.arange(lon_s, lon_e+1.e-6, dlon)
    lato1d = np.arange(lat_s, lat_e+1.e-6, dlat)
    nxo = len(lono1d)
    nyo = len(lato1d)
    lono, lato = np.meshgrid(lono1d, lato1d)

    X3dout = np.empty((nv3d, nz, nyo, nxo), dtype=X3d.dtype)
    X2dout = np.empty((nv2d, nyo, nxo), dtype=X2d.dtype)

    ri, rj = bmap(lono, lato)
    ri /= dx
    rj /= dy
    rij_interp = np.empty((nxo*nyo, 2), dtype=X3d.dtype)
    rij_interp[:,0] = rj.ravel()
    rij_interp[:,1] = ri.ravel()

    xic = np.arange(nx, dtype=X3d.dtype) + sio.bufsize
    xjc = np.arange(ny, dtype=X3d.dtype) + sio.bufsize

    missing = np.array(missing, dtype=X3d.dtype).ravel()[0]

    for iv in range(nv3d):
        for iz in range(nz):
            tmp2d = np.copy(X3d[iv,iz,:,:])
            tmp2d[tmp2d == missing] = np.nan
            X3dout[iv,iz,:,:] = interpn((xjc, xic), tmp2d, rij_interp, method='linear', bounds_error=False, fill_value=missing).reshape(nyo, nxo)

    for iv in range(nv2d):
        tmp2d = np.copy(X2d[iv,:,:])
        tmp2d[tmp2d == missing] = np.nan
        X2dout[iv,:,:] = interpn((xjc, xic), tmp2d, rij_interp, method='linear', bounds_error=False, fill_value=missing).reshape(nyo, nxo)

    return X3dout, X2dout, lon_s, lon_e, nxo, lat_s, lat_e, nyo


def create_ctlfile(sio, conf, nx, ny, nz, t, tint, tskip_a, nto_a, gradsfile, ctlfile):
    """
    """
    if sio.bufsize == 0:
        sliceobj = slice(None)
    else:
        sliceobj = slice(sio.bufsize, -sio.bufsize)
    lons = np.min(sio.lon[sliceobj])
    lone = np.max(sio.lon[sliceobj])
    lats = np.min(sio.lat[sliceobj])
    late = np.max(sio.lat[sliceobj])
    nxout = nx * 2 - 1
    nyout = ny * 2 - 1

    if conf['proj']['type'] == 'MER':
        nxout = nx
        nyout = ny
        merlat1 = sio.lat[sliceobj]
        merlat = ''
        for j in range(ny):
            merlat += "{0:12.6f}\n".format(merlat1[j,1])

    lonint = (lone - lons) / (nxout-1)
    latint = (late - lats) / (nyout-1)
    dx = sio.dimdef['coor_g']['x'][sio.bufsize+1] - sio.dimdef['coor_g']['x'][sio.bufsize]
    dy = sio.dimdef['coor_g']['y'][sio.bufsize+1] - sio.dimdef['coor_g']['y'][sio.bufsize]

    levs = ''
    if conf['vcoor'] == 'z' or conf['vcoor'] == 'o' or conf['ftype'] == 'restart_sprd':
        for ilev in range(nz):
            levs += "{0:12.6f}\n".format(sio.z[ilev])
    elif conf['vcoor'] == 'p':
        for ilev in range(nz):
            levs += "{0:12.6f}\n".format(conf['plevels'][ilev] / 100.)
    if conf['ftype'] == 'restart' or conf['ftype'] == 'restart_sprd':
        ts = t
    else:
        ts = sio.t[0]
        if len(sio.t) > 1:
            tint = (sio.t[1] - sio.t[0]) * tskip_a

    tint_min = int(round(tint.total_seconds() / 60))
    if tint_min == 0:
        tint_min = 1

    if conf['proj']['type'] == 'LC':
        if 'basepoint_x' in conf['proj'] and conf['proj']['basepoint_x'] is None:
            iref = 0.5*(nx+1)
        else:
            iref = conf['proj']['basepoint_x'] / float(dx) + 0.5
        if 'basepoint_y' in conf['proj'] and conf['proj']['basepoint_y'] is None:
            jref = 0.5*(ny+1)
        else:
            jref = conf['proj']['basepoint_y'] / float(dy) + 0.5
        pdef = 'pdef {isize:6d} {jsize:6d} lcc {latref:12.6f} {lonref:12.6f} {iref:.1f} {jref:.1f} {Struelat:12.6f} {Ntruelat:12.6f} {slon:12.6f} {dx:12.6f} {dy:12.6f}'.format(
               isize=nx, jsize=ny,
               latref=conf['proj']['basepoint_lat'], lonref=conf['proj']['basepoint_lon'],
               iref=iref, jref=jref,
               Struelat=conf['proj']['LC_lat1'], Ntruelat=conf['proj']['LC_lat2'],
               slon=conf['proj']['basepoint_lon'], dx=dx, dy=dy)
    elif conf['proj']['type'] == 'MER':
        if 'basepoint_x' in conf['proj'] and conf['proj']['basepoint_x'] is None:
            iref = 0.5*(nx+1)
        else:
            iref = conf['proj']['basepoint_x'] / float(dx) + 0.5
        if 'basepoint_y' in conf['proj'] and conf['proj']['basepoint_y'] is None:
            jref = 0.5*(ny+1)
        else:
            jref = conf['proj']['basepoint_y'] / float(dy) + 0.5
    else:
        raise ValueError('[Error] Unsupport map projection.')

    varstr = ''
    for ivar in conf['varout_3d']:
        varstr += "{varname:<12}{nz:6d} 99 {dscr:s}\n".format(varname=ivar, nz=nz, dscr=var_3d_name[ivar])
    for ivar in conf['varout_2d']:
        varstr += "{varname:<12}{nz:6d} 99 {dscr:s}\n".format(varname=ivar, nz=0, dscr=var_2d_name[ivar])

    if conf['proj']['type'] == 'LC':
      template = """dset ^{dset:s}
undef {undef:e}
xdef {nxout:6d} linear {lons:12.6f} {lonint:12.6f}
ydef {nyout:6d} linear {lats:12.6f} {latint:12.6f}
zdef {nz:6d} levels
{levs:s}tdef {nto:6d} linear {ts:s} {tint:d}mn
{pdef:s}
vars {nvar:d}
{varstr:s}endvars
"""
      context = {
      'dset':   os.path.relpath(gradsfile, os.path.dirname(ctlfile)),
      'undef':  conf['missing'],
      'nxout':  nxout,
      'lons':   lons,
      'lonint': lonint,
      'nyout':  nyout,
      'lats':   lats,
      'latint': latint,
      'nz':     nz,
      'levs':   levs,
      'nto':    nto_a,
      'ts':     ts.strftime('%H:%MZ%d%b%Y'),
      'tint':   tint_min,
      'pdef':   pdef,
      'nvar':   len(conf['varout_3d']) + len(conf['varout_2d']),
      'varstr': varstr
      }

    if conf['proj']['type'] == 'MER':
      template = """dset ^{dset:s}
undef {undef:e}
xdef {nxout:6d} linear {lons:12.6f} {lonint:12.6f}
ydef {nyout:6d} levels 
{merlat:s} 
zdef {nz:6d} levels
{levs:s}tdef {nto:6d} linear {ts:s} {tint:d}mn
vars {nvar:d}
{varstr:s}endvars
"""
      context = {
      'dset':   os.path.relpath(gradsfile, os.path.dirname(ctlfile)),
      'undef':  conf['missing'],
      'nxout':  nxout,
      'lons':   lons,
      'lonint': lonint,
      'nyout':  nyout,
      'lats':   lats,
      'latint': latint,
      'nz':     nz,
      'levs':   levs,
      'nto':    nto_a,
      'ts':     ts.strftime('%H:%MZ%d%b%Y'),
      'tint':   tint_min,
      'nvar':   len(conf['varout_3d']) + len(conf['varout_2d']),
      'varstr': varstr,
      'merlat': merlat
      }


    with open(ctlfile, 'w') as fc:
        fc.write(template.format(**context))


def create_ctlfile_ll(sio, conf, lons, lone, nx, lats, late, ny, nz, t, tint, tskip_a, nto_a, gradsfile, ctlfile):
    """
    """
    lonint = (lone - lons) / (nx-1)
    latint = (late - lats) / (ny-1)

    levs = ''
    if conf['vcoor'] == 'z' or conf['vcoor'] == 'o' or conf['ftype'] == 'restart_sprd':
        for ilev in range(nz):
            levs += "{0:12.6f}\n".format(sio.z[ilev])
    elif conf['vcoor'] == 'p':
        for ilev in range(nz):
            levs += "{0:12.6f}\n".format(conf['plevels'][ilev] / 100.)

    if conf['ftype'] == 'restart' or conf['ftype'] == 'restart_sprd':
        ts = t
    else:
        ts = sio.t[0]
        if len(sio.t) > 1:
            tint = (sio.t[1] - sio.t[0]) * tskip_a
    tint_min = int(round(tint.total_seconds() / 60))
    if tint_min == 0:
        tint_min = 1

    varstr = ''
    for ivar in conf['varout_3d']:
        varstr += "{varname:<12}{nz:6d} 99 {dscr:s}\n".format(varname=ivar, nz=nz, dscr=var_3d_name[ivar])
    for ivar in conf['varout_2d']:
        varstr += "{varname:<12}{nz:6d} 99 {dscr:s}\n".format(varname=ivar, nz=0, dscr=var_2d_name[ivar])

    template = """dset ^{dset:s}
undef {undef:e}
xdef {nx:6d} linear {lons:12.6f} {lonint:12.6f}
ydef {ny:6d} linear {lats:12.6f} {latint:12.6f}
zdef {nz:6d} levels
{levs:s}tdef {nto:6d} linear {ts:s} {tint:d}mn
vars {nvar:d}
{varstr:s}endvars
"""
    context = {
    'dset':   os.path.relpath(gradsfile, os.path.dirname(ctlfile)),
    'undef':  conf['missing'],
    'nx':     nx,
    'lons':   lons,
    'lonint': lonint,
    'ny':     ny,
    'lats':   lats,
    'latint': latint,
    'nz':     nz,
    'levs':   levs,
    'nto':    nto_a,
    'ts':     ts.strftime('%H:%MZ%d%b%Y'),
    'tint':   tint_min,
    'nvar':   len(conf['varout_3d']) + len(conf['varout_2d']),
    'varstr': varstr
    }

    with open(ctlfile, 'w') as fc:
        fc.write(template.format(**context))


def convert(basename, topo=None, t=dt.datetime(2000, 1, 1), tint=dt.timedelta(hours=6), 
            gradsfile='out.dat', ctlfile='auto', gradsfile_ll='out_latlon.dat', ctlfile_ll='auto',
            comm=None, commL=None, sim_read=1, **kwargs):
    """
    """
    # Determine if the mpi4py is used
    #------------
    if comm is None:
        nprocs = 1
        myrank = 0
        myrankmsg = ''
        nprocsL = 1
        myrankL = 0
    else:
        nprocs = comm.Get_size()
        myrank = comm.Get_rank()
        myrankmsg = '<< Rank {:6d} >> '.format(myrank)
        if commL is None:
            commL = comm
            print('<< My rank / total processes = {:d} / {:d} >>'.format(myrank, nprocs))
        commL.Barrier()
        nprocsL = commL.Get_size()
        myrankL = commL.Get_rank()

    # Initial settings
    #------------
    conf = rc(**kwargs)

    if conf['ftype'] in ('history', 'history_z'):
        bufsize = 0
    else:
        bufsize = 2
    if conf['ftype'] == 'restart_sprd':
        conf['varout_3d'] = [i for i in conf['varout_3d'] if i in var_3d_sprd]
        conf['varout_2d'] = [i for i in conf['varout_2d'] if i in var_2d_sprd]
    else:
        conf['varout_3d'] = [i for i in conf['varout_3d'] if i in var_3d]
        conf['varout_2d'] = [i for i in conf['varout_2d'] if i in var_2d]

    if sim_read <= 0 or sim_read >= nprocsL:
        sim_read = nprocsL

    # Initialize the ScaleIO object and get dimensions using the master process,
    # then broadcast the dimensions to other processes
    #------------
    if myrankL == 0:
        sio = ScaleIO(basename, cache=True, bufsize=bufsize, verbose=1)

        nx = sio.dimdef['len_g']['x'] - sio.bufsize * 2
        ny = sio.dimdef['len_g']['y'] - sio.bufsize * 2
        nz = sio.dimdef['len']['z'][0]
        if conf['vcoor'] == 'z' or conf['vcoor'] == 'o' or conf['ftype'] == 'restart_sprd':
            nzout = nz
        elif conf['vcoor'] == 'p':
            nzout = len(conf['plevels'])
        else:
            raise ValueError("vcoor = '{0:s}' is not supported. vcoor: {'z', 'p', 'o'}".format(conf['vcoor']))
        if sio.t is None:
            nt = 1
        else:
            nt = len(sio.t)
        bmap = set_bmap(sio, conf['proj'])

        print('--------------------')
        print('nx =', nx)
        print('ny =', ny)
        print('nz =', nz)
        print('nzout =', nzout)
        print('nt =', nt)
        print('--------------------')
    else:
        sio = None
        nx = 0
        ny = 0
        nz = 0
        nzout = 0
        nt = 0
        bmap = None

    if nprocsL > 1:
        nx = commL.bcast(nx, root=0)
        ny = commL.bcast(ny, root=0)
        nz = commL.bcast(nz, root=0)
        nzout = commL.bcast(nzout, root=0)
        nt = commL.bcast(nt, root=0)
        bmap = commL.bcast(bmap, root=0)

    # Determine the time frames handled by each process
    #------------
    its_a = conf['tstart']
    ite_a = conf['tend']
    tskip_a = conf['tskip']
    if ite_a == -1 or ite_a > nt:
        ite_a = nt
    nto_a = len(range(its_a, ite_a, tskip_a))

    its = its_a + tskip_a * myrankL
    ite = ite_a
    tskip = tskip_a * nprocsL
    nto = len(range(its, ite, tskip))

    # Generate the CTL file using the master process
    #------------
    if myrankL == 0:
        if gradsfile is not None and ctlfile is not None:
            print('Generate CTL file')
            if ctlfile is 'auto':
                ctlfile_ = '{0:s}.ctl'.format(gradsfile.rsplit('.', 1)[0])
            else:
                ctlfile_ = ctlfile
            create_ctlfile(sio, conf, nx, ny, nzout, t, tint, tskip_a, nto_a, gradsfile, ctlfile_)

        if gradsfile_ll is not None and ctlfile_ll is not None:
            lon = sio.readvar('lon')
            lat = sio.readvar('lat')
            lon_s = np.floor(np.min(lon) / conf['dlon']) * conf['dlon']
            lon_e = np.ceil(np.max(lon) / conf['dlon']) * conf['dlon']
            lat_s = np.floor(np.min(lat) / conf['dlat']) * conf['dlat']
            lat_e = np.ceil(np.max(lat) / conf['dlat']) * conf['dlat']
            nxll = int(np.rint((lon_e - lon_s) / conf['dlon'])) + 1
            nyll = int(np.rint((lat_e - lat_s) / conf['dlat'])) + 1

            print('Generate CTL file (lat/lon)')
            if ctlfile_ll is 'auto':
                ctlfile_ll_ = '{0:s}.ctl'.format(gradsfile_ll.rsplit('.', 1)[0])
            else:
                ctlfile_ll_ = ctlfile_ll
            create_ctlfile_ll(sio, conf, lon_s, lon_e, nxll, lat_s, lat_e, nyll,
                              nzout, t, tint, tskip_a, nto_a, gradsfile_ll, ctlfile_ll_)

    # Main computation
    #------------
    if nto > 0:
        var_necessary = []
        if conf['vcoor'] == 'z' and conf['ftype'] != 'restart_sprd':
            var_necessary += ['z']
            if conf['extrap']:
                var_necessary += ['p', 'tk', 'theta', 'rho', 'rhot', 'qv', 'qhydro']

        if conf['vcoor'] == 'p' and conf['ftype'] != 'restart_sprd':
            var_necessary += ['p']
            if conf['extrap']:
                var_necessary += ['z', 'tk', 'theta', 'rho', 'rhot', 'qv', 'qhydro']
#            if 'qhydro' in conf['varout_3d']:
#                var_necessary += ['qc', 'qr', 'qi', 'qs', 'qg']
            if 'dbz' in conf['varout_3d'] or 'max_dbz' in conf['varout_2d']:
                var_necessary += ['rho', 'qr', 'qs', 'qg']
        if 'u' in conf['varout_3d'] or 'v' in conf['varout_3d']:
            var_necessary += ['u', 'v']
        if 'u10' in conf['varout_2d'] or 'v10' in conf['varout_2d']:
            var_necessary += ['u10', 'v10']
        if conf['ftype'] == 'restart':
            if 'psfc' in conf['varout_2d']:
                var_necessary += ['rho', 'p', 'z']

        dummy = np.zeros(1)

        ito = 0
        for it in range(its, ite, tskip):
            it_a = its_a + tskip * ito + tskip_a * myrankL
            ito_a = nprocsL * ito + myrankL

            # read data in 'dryrun' mode, one process follows the previous process sequentially
            ######
            if commL is not None and ito_a >= sim_read:
                srank = myrankL - sim_read
                if srank < 0:
                    srank += nprocsL
#                print(myrankmsg, 'Recv', 10+ito_a-sim_read)
                commL.Recv(dummy, source=srank, tag=10+ito_a-sim_read)
            ######
            if sio is None:
                sio = ScaleIO(basename, cache=True, bufsize=bufsize, verbose=1)
            if conf['ftype'] != 'restart_sprd':
                if topo is None:
                    if conf['ftype'] == 'restart':
                        topo = sio.readvar('TOPO')
                    elif conf['ftype'] == 'history':
                        topo = sio.readvar('topo', t=0)
                elif type(topo) is str :
                    sio_topo = ScaleIO(topo, cache=True, bufsize=2, verbose=1)
                    topo = sio_topo.readvar('TOPO')
                    del sio_topo
            X, t0_ext = convert_readvar(sio, bmap, topo, conf, var_necessary, it, tskip_a, dryrun=True)
            ######
            if commL is not None and ito_a + sim_read < nto_a:
                drank = myrankL + sim_read
                if drank >= nprocsL:
                    drank -= nprocsL
#                print(myrankmsg, 'Send', 10+ito_a)
                commL.Send(dummy, dest=drank, tag=10+ito_a)
            ######

            # do real variable transform computation
            X, t0_ext = convert_readvar(sio, bmap, topo, conf, var_necessary, it, tskip_a, myrankmsg=myrankmsg)

            # vertical interpolation and extrapolation
            X3d, X2d, nv3d, nv2d = convert_vintp(sio, conf, X, t0_ext, nx, ny, nzout, it, myrankmsg=myrankmsg)

            # horizontal interpolation
            if gradsfile_ll is not None:
                print(myrankmsg, 'Convert to lat/lon grid [to = {:d}]'.format(ito_a+1))
                X3dll, X2dll, lons, lone, nxout, lats, late, nyout = \
                    convert_hintp(sio, bmap, conf, X3d, X2d, conf['dlon'], conf['dlat'], conf['missing'])

            sio.freecache()

            if gradsfile is not None or gradsfile_ll is not None:
                # write grads data, one process follows the previous process sequentially
                ######
                if commL is not None and ito_a >= sim_read:
                    srank = myrankL - sim_read
                    if srank < 0:
                        srank += nprocsL
#                    print(myrankmsg, 'Recv', 10+nto_a+ito_a-sim_read)
                    commL.Recv(dummy, source=srank, tag=10+nto_a+ito_a-sim_read)
                ######
                if gradsfile is not None:
                    if ito_a == 0:
                        f = open(gradsfile, 'wb')
                    else:
                        while not os.path.exists(gradsfile):
                            time.sleep(1)
                        f = open(gradsfile, 'r+b')
                    for iv in range(nv3d):
                        print(myrankmsg, 'Write 3D variable: {:s} [to = {:d}]'.format(conf['varout_3d'][iv], ito_a+1))
                        gradsio.writegrads(f, X3d[iv], iv+1, nv3d=nv3d, nv2d=nv2d, t=ito_a+1, nx=nx, ny=ny, nz=nzout, nt=nt)
                    for iv in range(nv2d):
                        print(myrankmsg, 'Write 2D variable: {:s} [to = {:d}]'.format(conf['varout_2d'][iv], ito_a+1))
                        gradsio.writegrads(f, X2d[iv], nv3d+iv+1, nv3d=nv3d, nv2d=nv2d, t=ito_a+1, nx=nx, ny=ny, nz=nzout, nt=nt)
                    f.close()

                if gradsfile_ll is not None:
                    if ito_a == 0:
                        f2 = open(gradsfile_ll, 'wb')
                    else:
                        while not os.path.exists(gradsfile_ll):
                            time.sleep(1)
                        f2 = open(gradsfile_ll, 'r+b')
                    for iv in range(nv3d):
                        print(myrankmsg, 'Write 3D variable (lat/lon): {:s} [to = {:d}]'.format(conf['varout_3d'][iv], ito_a+1))
                        gradsio.writegrads(f2, X3dll[iv], iv+1, nv3d=nv3d, nv2d=nv2d, t=ito_a+1, nx=nxout, ny=nyout, nz=nzout, nt=nt)
                    for iv in range(nv2d):
                        print(myrankmsg, 'Write 2D variable (lat/lon): {:s} [to = {:d}]'.format(conf['varout_2d'][iv], ito_a+1))
                        gradsio.writegrads(f2, X2dll[iv], nv3d+iv+1, nv3d=nv3d, nv2d=nv2d, t=ito_a+1, nx=nxout, ny=nyout, nz=nzout, nt=nt)
                    f2.close()
                ######
                if commL is not None and ito_a + sim_read < nto_a:
                    drank = myrankL + sim_read
                    if drank >= nprocsL:
                        drank -= nprocsL
#                    print(myrankmsg, 'Send', 10+nto_a+ito_a)
                    commL.Send(dummy, dest=drank, tag=10+nto_a+ito_a)
                ######

            X.clear()

            ito += 1

    if sio is not None:
        del sio
