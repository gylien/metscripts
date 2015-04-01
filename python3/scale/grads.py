import numpy as np
#import numpy.ma as ma
import datetime as dt
from .io import ScaleIO
from .calc import *
import gradsio


__all__ = ['convert']


missing = -9.99e33
vcoor = 'p'

var_3d = {
'u': 'u-wind (m/s)',
'v': 'v-wind (m/s)',
'w': 'w-wind (m/s)',
'p': 'Pressure (Pa)',
't': 'Temperature (K)',
'theta': 'Potential temperature (K)',
'rho': 'Density (kg/m^3)',
'qv': 'Water vapor mixing ratio (kg/kg)',
'qc': 'Cloud water mixing ratio (kg/kg)',
'qr': 'Rain water mixing ratio (kg/kg)',
'qi': 'Cloud ice mixing ratio (kg/kg)',
'qs': 'Snow mixing ratio (kg/kg)',
'qg': 'Graupel mixing ratio (kg/kg)',
'qhydro': 'Mixing ratio of all hydrometers (kg/kg)',
'dbz': 'Radar reflectivity (dBZ)'
}
var_2d = {
'topo': 'Topography height (m)',
'max_dbz': 'Maximum radar reflectivity (dBZ)'
}

varout_3d = ['u', 'v', 'w', 'p', 't', 'theta', 'rho', 'qv', 'qc', 'qr', 'qi', 'qs', 'qg', 'qhydro', 'dbz']
#
#
varout_2d = ['topo', 'max_dbz']
#
# slp
plevels = [100000., 92500., 85000., 70000., 50000., 30000., 20000., 10000., 5000., 2000., 1000.]

proj = {
'type': 'LC',
'basepoint_lon': 135.220404,
'basepoint_lat': 34.653396,
'LC_lat1': 30.0,
'LC_lat2': 40.0,
}


def convert(basename, basename_topo=None, ftype='restart', gradsfile='out', t=dt.datetime(2000, 1, 1)):
    """
    """
    sio = ScaleIO(basename)
    nx = sio.dimdef['len_g']['x']
    ny = sio.dimdef['len_g']['y']
    if ftype == 'restart':
        nx -= 4
        ny -= 4
    nz = sio.dimdef['len']['z'][0]
    if vcoor == 'z':
        nzout = nz
    elif vcoor == 'p':
        nzout = len(plevels)
    else:
        raise ValueError("vcoor = '{0:s}' is not supported. vcoor: {'z', 'p'}".format(vcoor))
    if sio.t is None:
        nt = 1
    else:
        nt = len(sio.t)
    print('nx =', nx)
    print('ny =', ny)
    print('nz =', nz, ', nzout =', nzout)
    print('nt =', nt)
    print()

    var = {}
    if basename_topo is None:
        var['topo'] = sio.readvar('TOPO')[2:-2,2:-2]
    else:
        sio_topo = ScaleIO(basename_topo)
        var['topo'] = sio_topo.readvar('TOPO')[2:-2,2:-2]
#    var['topo'] = np.tile(var['topo'], (nt, 1, 1))

    f = open(gradsfile + '.dat', 'wb')

    for it in range(nt):
        if ftype == 'restart':
            var['rho'] = sio.readvar('DENS', t=it)[:,2:-2,2:-2]
            var['u'], var['v'], var['w'] = calc_uvw(sio, rho=var['rho'], first_grd=True, t=it)

            var['qv'] = sio.readvar('QV', t=it)[:,2:-2,2:-2]
            var['qc'] = sio.readvar('QC', t=it)[:,2:-2,2:-2]
            var['qr'] = sio.readvar('QR', t=it)[:,2:-2,2:-2]
            var['qi'] = sio.readvar('QI', t=it)[:,2:-2,2:-2]
            var['qs'] = sio.readvar('QS', t=it)[:,2:-2,2:-2]
            var['qg'] = sio.readvar('QG', t=it)[:,2:-2,2:-2]
            var['qhydro'] = var['qc'] + var['qr'] + var['qi'] + var['qs'] + var['qg']
            var['p'], var['t'], var['theta'] = calc_pt(sio, rho=var['rho'], qv=var['qv'], qhydro=var['qhydro'], tout=True, thetaout=True, t=it)

            if 'dbz' in varout_3d or 'max_dbz' in varout_2d:
                var['dbz'], var['max_dbz'] = calc_ref(sio, rho=var['rho'], qr=var['qr'], qs=var['qs'], qg=var['qg'], t=it)

            if vcoor == 'z':
                height, height_h = calc_height(sio, topo=var['topo'])
                for ivar in varout_3d:
                    print('Vertical interpolation at Z-coordinate: ', ivar)
                    var[ivar] = interp_z(sio, var[ivar], height=height, t=it)
        elif ftype == 'history':
            print("Not finished...")
        else:
            raise ValueError("ftype = '{0:s}' is not supported. ftype: {'restart', 'history'}".format(ftype))

        if vcoor == 'p':
            for ivar in varout_3d:
                print('Vertical interpolation at P-coordinate: ', ivar)
                interp_p(sio, var[ivar], plevels, p=var['p'], t=it)


        for iv, ivar in enumerate(varout_3d):
            gradsio.writegrads(f, var[ivar].filled(fill_value=missing), iv+1,
                               nv3d=len(varout_3d), nv2d=len(varout_2d), t=it+1, nx=nx, ny=ny, nz=nz, nt=nt)
        for iv, ivar in enumerate(varout_2d):
            gradsio.writegrads(f, var[ivar].filled(fill_value=missing), len(varout_3d)+iv+1,
                               nv3d=len(varout_3d), nv2d=len(varout_2d), t=it+1, nx=nx, ny=ny, nz=nz, nt=nt)

    f.close()

    if ftype == 'restart':
        lons = np.min(sio.lon[2:-2,2:-2])
        lone = np.max(sio.lon[2:-2,2:-2])
        lats = np.min(sio.lat[2:-2,2:-2])
        late = np.max(sio.lat[2:-2,2:-2])
    else:
        lons = np.min(sio.lon)
        lone = np.max(sio.lon)
        lats = np.min(sio.lat)
        late = np.max(sio.lat)
    nxout = nx * 2 - 1
    nyout = ny * 2 - 1
    lonint = (lone - lons) / (nxout-1)
    latint = (late - lats) / (nyout-1)
    if ftype == 'restart':
        dx = sio.dimdef['coor_g']['x'][3] - sio.dimdef['coor_g']['x'][2]
        dy = sio.dimdef['coor_g']['y'][3] - sio.dimdef['coor_g']['y'][2]
    else:
        dx = sio.dimdef['coor_g']['x'][1] - sio.dimdef['coor_g']['x'][0]
        dy = sio.dimdef['coor_g']['y'][0] - sio.dimdef['coor_g']['y'][0]

    levs = ''
    if vcoor == 'z':
        for ilev in range(len(sio.z)):
            levs += "{0:12.6f}\n".format(sio.z[ilev])
    elif vcoor == 'p':
        for ilev in range(len(plevels)):
            levs += "{0:12.6f}\n".format(plevels[ilev] / 100.)

    if ftype == 'restart':
        ts = t
        tint = dt.timedelta(hours=6)
    else:
        ts = sio.t[0]
        if len(sio.t) > 1:
            tint = sio.t[1] - sio.t[0]
        else:
            tint = dt.timedelta(hours=6)

    pdef = 'pdef {isize:6d} {jsize:6d} lccr {latref:12.6f} {lonref:12.6f} {iref:.1f} {jref:.1f} {Struelat:12.6f} {Ntruelat:12.6f} {slon:12.6f} {dx:12.6f} {dy:12.6f}'.format(
           isize=nx, jsize=ny,
           latref=proj['basepoint_lat'], lonref=proj['basepoint_lon'],
           iref=0.5*(nx+1), jref=0.5*(ny+1),
           Struelat=proj['LC_lat1'], Ntruelat=proj['LC_lat2'],
           slon=proj['basepoint_lon'], dx=dx, dy=dy)

    var = ''
    for ivar in varout_3d:
        var += "{varname:<12}{nz:6d} 99 {dscr:s}\n".format(varname=ivar, nz=nzout, dscr=var_3d[ivar])
    for ivar in varout_2d:
        var += "{varname:<12}{nz:6d} 99 {dscr:s}\n".format(varname=ivar, nz=0, dscr=var_2d[ivar])

    template = """dset ^{dset:s}.dat
undef {undef:e}
xdef {nxout:6d} linear {lons:12.6f} {lonint:12.6f}
ydef {nyout:6d} linear {lats:12.6f} {latint:12.6f}
zdef {nz:6d} levels
{levs:s}tdef {nt:6d} linear {ts:s} {tint:s}
{pdef:s}
vars {nvar:d}
{var:s}endvars
"""
    context = {
    'dset':   gradsfile,
    'undef':  missing,
    'nxout':  nxout,
    'lons':   lons,
    'lonint': lonint,
    'nyout':  nyout,
    'lats':   lats,
    'latint': latint,
    'nz':     nzout,
    'levs':   levs,
    'nt':     nt,
    'ts':     ts.strftime('%H:%MZ%d%b%Y'),
    'tint':   '{0:d}mn'.format(int(round(tint.total_seconds() / 60)))	,
    'pdef':   pdef,
    'nvar':   len(varout_3d) + len(varout_2d),
    'var':    var,
    }

    with open(gradsfile + '.ctl', 'w') as fc:
        fc.write(template.format(**context))



