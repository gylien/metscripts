import numpy as np
import numpy.ma as ma
import datetime as dt
import os
from .io import ScaleIO
from .calc import *
import gradsio


__all__ = ['convert']


config = {
'ftype': 'restart',
'missing': -9.99e33,
'vcoor': 'p',
'plevels': [100000., 92500., 85000., 70000., 50000., 30000., 20000., 10000., 5000., 2000., 1000.],
'varout_3d': ['u', 'v', 'w', 'p', 'temp', 'theta', 'rho', 'momx', 'momy', 'momz', 'rhot', 'z', 'qv', 'qc', 'qr', 'qi', 'qs', 'qg', 'qhydro', 'dbz'],
'varout_2d': ['topo', 'max_dbz'],
'proj': {'type': 'LC',
         'basepoint_lon': 136.,
         'basepoint_lat': 36.,
         'LC_lat1': 30.0,
         'LC_lat2': 42.0,}
}

#proj = {
#'type': 'LC',
#'basepoint_lon': 135.220404,
#'basepoint_lat': 34.653396,
#'LC_lat1': 30.0,
#'LC_lat2': 40.0,
#}
#proj = {
#'type': 'LC',
#'basepoint_lon': 135.52,
#'basepoint_lat': 34.82,
#'LC_lat1': 34.52,
#'LC_lat2': 35.12,
#}

var_3d = {
'u': 'u-wind (m/s)',
'v': 'v-wind (m/s)',
'w': 'w-wind (m/s)',
'p': 'Pressure (Pa)',
'temp': 'Temperature (K)',
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


def rc(**kwargs):
    for key, value in list(kwargs.items()):
        if key in config:
            if value != '-':
                config[key] = value
        else:
            raise KeyError("'{0:s}' is not a configuration key.".format(key))


def convert(basename, topo=None, gradsfile='out.dat', ctlfile='auto', t=dt.datetime(2000, 1, 1), **kwargs):
    """
    """
    rc(**kwargs)

    sio = ScaleIO(basename)
    nx = sio.dimdef['len_g']['x']
    ny = sio.dimdef['len_g']['y']
    if config['ftype'] == 'restart' or config['ftype'] == 'restart_sprd':
        nx -= 4
        ny -= 4
    nz = sio.dimdef['len']['z'][0]
    if config['vcoor'] == 'z' or config['vcoor'] == 'o' or config['ftype'] == 'restart_sprd':
        nzout = nz
    elif config['vcoor'] == 'p':
        nzout = len(config['plevels'])
    else:
        raise ValueError("vcoor = '{0:s}' is not supported. vcoor: {'z', 'p', 'o'}".format(config['vcoor']))
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
    if config['ftype'] != 'restart_sprd':
        if topo is None:
            print('Read variable: TOPO')
            if config['ftype'] == 'restart':
                var['topo'] = sio.readvar('TOPO')[2:-2,2:-2]
            elif config['ftype'] == 'history':
                var['topo'] = sio.readvar('topo')
        elif type(topo) is str :
            sio_topo = ScaleIO(topo)
            print('Read variable: TOPO')
            var['topo'] = sio_topo.readvar('TOPO')[2:-2,2:-2]
        else:
            var['topo'] = topo

    f = open(gradsfile, 'wb')

    for it in range(nt):
        if config['ftype'] == 'restart':
            for ivar, ivarf in ('rho', 'DENS'), ('momz', 'MOMZ'), ('rhot', 'RHOT'), \
                               ('qv', 'QV'), ('qc', 'QC'), ('qr', 'QR'), ('qi', 'QI'), ('qs', 'QS'), ('qg', 'QG'):
                print('Read variable: ' + ivarf)
                var[ivar] = sio.readvar(ivarf, t=it)[:,2:-2,2:-2]
            print('Read variable: MOMX')
            var['momx'] = sio.readvar('MOMX', t=it)[:,2:-2,1:-2]
            print('Read variable: MOMY')
            var['momy'] = sio.readvar('MOMY', t=it)[:,1:-2,2:-2]

            print('Calculate: u, v, w')
            var['u'], var['v'], var['w'], var['momx'], var['momy'], var['momz'] = \
                calc_uvw(sio, rho=var['rho'], momx=var['momx'], momy=var['momy'], momz=var['momz'], first_grd=True, t=it)
            print('Calculate: qhydro')
            var['qhydro'] = var['qc'] + var['qr'] + var['qi'] + var['qs'] + var['qg']
            print('Calculate: p, t, theta')
            var['p'], var['temp'], var['theta'] = calc_pt(sio, rho=var['rho'], rhot=var['rhot'], qv=var['qv'], qhydro=var['qhydro'], tout=True, thetaout=True, t=it)
            if 'dbz' in config['varout_3d'] or 'max_dbz' in config['varout_2d']:
                print('Calculate: dbz, max_dbz')
                var['dbz'], var['max_dbz'] = calc_ref(sio, rho=var['rho'], qr=var['qr'], qs=var['qs'], qg=var['qg'], t=it)
            if 'z' in config['varout_3d'] or config['vcoor'] == 'z':
                print('Calculate: z')
                var['z'], height_h = calc_height(sio, topo=var['topo'])

            if config['vcoor'] == 'z':
                for ivar in config['varout_3d']:
                    print('Vertical interpolation at Z-coordinate: ', ivar)
                    var_itp[ivar] = interp_z(sio, var[ivar], height=var['z'], t=it)

        elif config['ftype'] == 'restart_sprd':
            for ivar, ivarf in ('u', 'DENS'), ('v', 'MOMX'), ('w', 'MOMY'), ('temp', 'MOMZ'), ('p', 'RHOT'), \
                               ('qv', 'QV'), ('qc', 'QC'), ('qr', 'QR'), ('qi', 'QI'), ('qs', 'QS'), ('qg', 'QG'):
                if ivar in config['varout_3d']:
                    print('Read variable: ' + ivarf)
                    var[ivar] = sio.readvar(ivarf, t=it)[:,2:-2,2:-2]

            print('Destagger: u, v, w')
            if 'u' in config['varout_3d']:
                var['u'] = calc_destagger(var['u'], axis=2, first_grd=False)
            if 'v' in config['varout_3d']:
                var['v'] = calc_destagger(var['v'], axis=1, first_grd=False)
            if 'w' in config['varout_3d']:
                var['w'] = calc_destagger(var['w'], axis=0, first_grd=False)
                var['w'][0,:,:] = 0.

        elif config['ftype'] == 'history':
            for ivar, ivarf in ('rho', 'DENS'), ('momx', 'MOMX'), ('momy', 'MOMY'), ('momz', 'MOMZ'), ('rhot', 'RHOT'), \
                               ('qv', 'QV'), ('qc', 'QC'), ('qr', 'QR'), ('qi', 'QI'), ('qs', 'QS'), ('qg', 'QG'), \
                               ('u', 'U'), ('v', 'V'), ('w', 'W'), ('temp', 'T'), ('p', 'PRES'), ('theta', 'PT'):
                if ivar in config['varout_3d']:
                    print('Read variable: ' + ivarf)
                    var[ivar] = sio.readvar(ivarf, t=it)

            if 'qhydro' in config['varout_3d']:
                print('Calculate: qhydro')
                var['qhydro'] = var['qc'] + var['qr'] + var['qi'] + var['qs'] + var['qg']
            if 'dbz' in config['varout_3d'] or 'max_dbz' in config['varout_2d']:
                print('Calculate: dbz, max_dbz')
                var['dbz'], var['max_dbz'] = calc_ref(sio, rho=var['rho'], qr=var['qr'], qs=var['qs'], qg=var['qg'], t=it)
            if 'z' in config['varout_3d']:
                var['z'] = np.empty((len(sio.z), ny, nx), dtype=sio.z.dtype)
                for ilev in range(len(sio.z)):
                    var['z'][ilev] = sio.z[ilev]

        else:
            raise ValueError("ftype = '{0:s}' is not supported. ftype: {'restart', 'restart_sprd', 'history'}".format(config['ftype']))

        if config['vcoor'] == 'p':
            for ivar in config['varout_3d']:
                if ivar == 'p':
                    varshape = list(var[ivar].shape)
                    varshape[0] = len(config['plevels'])
                    var_itp[ivar] = np.empty(varshape, dtype=var[ivar].dtype)
                    for ilev in range(len(config['plevels'])):
                        var_itp[ivar][ilev] = config['plevels'][ilev]
                else:
                    print('Vertical interpolation at P-coordinate: ', ivar)
                    var_itp[ivar] = interp_p(sio, var[ivar], config['plevels'], p=var['p'], t=it)


        varout_3d_final = []
        varout_2d_final = []
        if config['vcoor'] == 'o' or (config['vcoor'] == 'z' and config['ftype'] == 'history') or config['ftype'] == 'restart_sprd':
            var_out = var
        else:
            var_out = var_itp

        for ivar in config['varout_3d']:
            if ivar in var_out:
                varout_3d_final.append(ivar)
        for ivar in config['varout_2d']:
            if ivar in var:  # always look for 'var' (instead of 'var_out') for 2D variables
                varout_2d_final.append(ivar)

        for iv, ivar in enumerate(varout_3d_final):
            print('Write 3D variable: {:s}'.format(ivar))
            if type(var_out[ivar]) == ma.MaskedArray:
                varf = var_out[ivar].filled(fill_value=config['missing'])
            else:
                varf = var_out[ivar]

            gradsio.writegrads(f, varf, iv+1,
                               nv3d=len(varout_3d_final), nv2d=len(varout_2d_final), t=it+1, nx=nx, ny=ny, nz=nzout, nt=nt)

        for iv, ivar in enumerate(varout_2d_final):
            print('Write 2D variable: {:s}'.format(ivar))
            if type(var[ivar]) == ma.MaskedArray:
                varf = var[ivar].filled(fill_value=config['missing'])  # always look for 'var' (instead of 'var_out') for 2D variables
            else:
                varf = var[ivar]
            gradsio.writegrads(f, varf, len(varout_3d_final)+iv+1,
                               nv3d=len(varout_3d_final), nv2d=len(varout_2d_final), t=it+1, nx=nx, ny=ny, nz=nzout, nt=nt)

    f.close()


    if ctlfile is not None:
        print('Generate CTL file')

        if config['ftype'] == 'restart' or config['ftype'] == 'restart_sprd':
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
        if config['ftype'] == 'restart' or config['ftype'] == 'restart_sprd':
            dx = sio.dimdef['coor_g']['x'][3] - sio.dimdef['coor_g']['x'][2]
            dy = sio.dimdef['coor_g']['y'][3] - sio.dimdef['coor_g']['y'][2]
        else:
            dx = sio.dimdef['coor_g']['x'][1] - sio.dimdef['coor_g']['x'][0]
            dy = sio.dimdef['coor_g']['y'][1] - sio.dimdef['coor_g']['y'][0]

        levs = ''
        if config['vcoor'] == 'z' or config['vcoor'] == 'o' or config['ftype'] == 'restart_sprd':
            for ilev in range(len(sio.z)):
                levs += "{0:12.6f}\n".format(sio.z[ilev])
        elif config['vcoor'] == 'p':
            for ilev in range(len(config['plevels'])):
                levs += "{0:12.6f}\n".format(config['plevels'][ilev] / 100.)

        if config['ftype'] == 'restart' or config['ftype'] == 'restart_sprd':
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
               latref=config['proj']['basepoint_lat'], lonref=config['proj']['basepoint_lon'],
               iref=0.5*(nx+1), jref=0.5*(ny+1),
               Struelat=config['proj']['LC_lat1'], Ntruelat=config['proj']['LC_lat2'],
               slon=config['proj']['basepoint_lon'], dx=dx, dy=dy)

        var = ''
        for ivar in varout_3d_final:
            var += "{varname:<12}{nz:6d} 99 {dscr:s}\n".format(varname=ivar, nz=nzout, dscr=var_3d[ivar])
        for ivar in varout_2d_final:
            var += "{varname:<12}{nz:6d} 99 {dscr:s}\n".format(varname=ivar, nz=0, dscr=var_2d[ivar])

        if ctlfile is 'auto':
            ctlfile0 = '{0:s}.ctl'.format(gradsfile.rsplit('.', 1)[0])
        else:
            ctlfile0 = ctlfile

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
        'dset':   os.path.relpath(gradsfile, os.path.dirname(ctlfile0)),
        'undef':  config['missing'],
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
        'nvar':   len(varout_3d_final) + len(varout_2d_final),
        'var':    var,
        }

        with open(ctlfile0, 'w') as fc:
            fc.write(template.format(**context))


