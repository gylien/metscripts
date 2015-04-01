import numpy as np
import numpy.ma as ma
import datetime as dt
from .io import ScaleIO
from .calc import *
import gradsio


__all__ = ['convert']


missing = -9.99e33
vcoor = 'o'

var_3d = {
'u': 'u-wind (m/s)',
'v': 'v-wind (m/s)',
'w': 'w-wind (m/s)',
'p': 'Pressure (Pa)',
't': 'Temperature (K)',
'theta': 'Potential temperature (K)',
'rho': 'Density (kg/m^3)',
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
'dbz': 'Radar reflectivity (dBZ)'
}
var_2d = {
'topo': 'Topography height (m)',
'max_dbz': 'Maximum radar reflectivity (dBZ)',
'slp': 'Sea level pressure (Pa)' ##### not yet implemented!!!
}

varout_3d = ['u', 'v', 'w', 'p', 't', 'theta', 'rho', 'momx', 'momy', 'momz', 'rhot', 'z', 'qv', 'qc', 'qr', 'qi', 'qs', 'qg', 'qhydro', 'dbz']
varout_2d = ['topo', 'max_dbz']

plevels = [100000., 92500., 85000., 70000., 50000., 30000., 20000., 10000., 5000., 2000., 1000.]

proj = {
'type': 'LC',
'basepoint_lon': 135.220404,
'basepoint_lat': 34.653396,
'LC_lat1': 30.0,
'LC_lat2': 40.0,
}
#proj = {
#'type': 'LC',
#'basepoint_lon': 135.52,
#'basepoint_lat': 34.82,
#'LC_lat1': 34.52,
#'LC_lat2': 35.12,
#}


def convert(basename, topo=None, ftype='restart', gradsfile='out.dat', ctlfile='auto', t=dt.datetime(2000, 1, 1)):
    """
    """
    sio = ScaleIO(basename)
    nx = sio.dimdef['len_g']['x']
    ny = sio.dimdef['len_g']['y']
    if ftype == 'restart':
        nx -= 4
        ny -= 4
    nz = sio.dimdef['len']['z'][0]
    if vcoor == 'z' or vcoor == 'o':
        nzout = nz
    elif vcoor == 'p':
        nzout = len(plevels)
    else:
        raise ValueError("vcoor = '{0:s}' is not supported. vcoor: {'z', 'p', 'o'}".format(vcoor))
    if sio.t is None:
        nt = 1
    else:
        nt = len(sio.t)
    print('--------------------')
    print('nx =', nx)
    print('ny =', ny)
    print('nz =', nz)
    print('nzout =', nzout)
    print('nt =', nt)
    print('--------------------')

    var = {}
    var_itp = {}
    if topo is None:
        print('Read variable: TOPO')
        var['topo'] = sio.readvar('TOPO')[2:-2,2:-2]
    elif type(topo) is str :
        sio_topo = ScaleIO(topo)
        print('Read variable: TOPO')
        var['topo'] = sio_topo.readvar('TOPO')[2:-2,2:-2]
    else:
        var['topo'] = topo

    f = open(gradsfile, 'wb')

    for it in range(nt):
        if ftype == 'restart':
            print('Read variable: DENS')
            var['rho'] = sio.readvar('DENS', t=it)[:,2:-2,2:-2]
            print('Read variable: MOMX')
            var['momx'] = sio.readvar('MOMX', t=it)[:,2:-2,1:-2]
            print('Read variable: MOMY')
            var['momy'] = sio.readvar('MOMY', t=it)[:,1:-2,2:-2]
            print('Read variable: MOMZ')
            var['momz'] = sio.readvar('MOMZ', t=it)[:,2:-2,2:-2]
            print('Read variable: RHOT')
            var['rhot'] = sio.readvar('RHOT', t=it)[:,2:-2,2:-2]
            print('Read variable: QV')
            var['qv'] = sio.readvar('QV', t=it)[:,2:-2,2:-2]
            print('Read variable: QC')
            var['qc'] = sio.readvar('QC', t=it)[:,2:-2,2:-2]
            print('Read variable: QR')
            var['qr'] = sio.readvar('QR', t=it)[:,2:-2,2:-2]
            print('Read variable: QI')
            var['qi'] = sio.readvar('QI', t=it)[:,2:-2,2:-2]
            print('Read variable: QS')
            var['qs'] = sio.readvar('QS', t=it)[:,2:-2,2:-2]
            print('Read variable: QG')
            var['qg'] = sio.readvar('QG', t=it)[:,2:-2,2:-2]

            print('Calculate: u, v, w')
            var['u'], var['v'], var['w'], var['momx'], var['momy'], var['momz'] = \
                calc_uvw(sio, rho=var['rho'], momx=var['momx'], momy=var['momy'], momz=var['momz'], first_grd=True, t=it)
            print('Calculate: qhydro')
            var['qhydro'] = var['qc'] + var['qr'] + var['qi'] + var['qs'] + var['qg']
            print('Calculate: p, t, theta')
            var['p'], var['t'], var['theta'] = calc_pt(sio, rho=var['rho'], rhot=var['rhot'], qv=var['qv'], qhydro=var['qhydro'], tout=True, thetaout=True, t=it)
            if 'dbz' in varout_3d or 'max_dbz' in varout_2d:
                print('Calculate: dbz, max_dbz')
                var['dbz'], var['max_dbz'] = calc_ref(sio, rho=var['rho'], qr=var['qr'], qs=var['qs'], qg=var['qg'], t=it)
            if 'z' in varout_3d or vcoor == 'z':
                print('Calculate: z')
                var['z'], height_h = calc_height(sio, topo=var['topo'])

            if vcoor == 'z':
                for ivar in varout_3d:
                    print('Vertical interpolation at Z-coordinate: ', ivar)
                    var_itp[ivar] = interp_z(sio, var[ivar], height=var['z'], t=it)
        elif ftype == 'history':
            print("Not finished...")
        else:
            raise ValueError("ftype = '{0:s}' is not supported. ftype: {'restart', 'history'}".format(ftype))

        if vcoor == 'p':
            for ivar in varout_3d:
                if ivar == 'p':
                    varshape = list(var[ivar].shape)
                    varshape[0] = len(plevels)
                    var_itp[ivar] = np.empty(varshape, dtype=var[ivar].dtype)
                    for ilev in range(len(plevels)):
                        var_itp[ivar][ilev] = plevels[ilev]
                else:
                    print('Vertical interpolation at P-coordinate: ', ivar)
                    var_itp[ivar] = interp_p(sio, var[ivar], plevels, p=var['p'], t=it)


        for iv, ivar in enumerate(varout_3d):
            print('Write 3D variable: {:s}'.format(ivar))
            if vcoor == 'o':
                if type(var[ivar]) == ma.MaskedArray:
                    varf = var[ivar].filled(fill_value=missing)
                else:
                    varf = var[ivar]
            else:
                if type(var_itp[ivar]) == ma.MaskedArray:
                    varf = var_itp[ivar].filled(fill_value=missing)
                else:
                    varf = var_itp[ivar]
            gradsio.writegrads(f, varf, iv+1,
                               nv3d=len(varout_3d), nv2d=len(varout_2d), t=it+1, nx=nx, ny=ny, nz=nzout, nt=nt)
        for iv, ivar in enumerate(varout_2d):
            print('Write 2D variable: {:s}'.format(ivar))
            if type(var[ivar]) == ma.MaskedArray:
                varf = var[ivar].filled(fill_value=missing)
            else:
                varf = var[ivar]
            gradsio.writegrads(f, varf, len(varout_3d)+iv+1,
                               nv3d=len(varout_3d), nv2d=len(varout_2d), t=it+1, nx=nx, ny=ny, nz=nzout, nt=nt)

    f.close()


    if ctlfile is not None:
        print('Generate CTL file')

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
        if vcoor == 'z' or vcoor == 'o':
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

        template = """dset ^{dset:s}
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

        if ctlfile is 'auto':
            ctlfile0 = '{0:s}.ctl'.format(gradsfile.rsplit('.', 1)[0])
        else:
            ctlfile0 = ctlfile
        with open(ctlfile0, 'w') as fc:
            fc.write(template.format(**context))


