import numpy as np
import numpy.ma as ma
import datetime as dt
import os
from mpl_toolkits.basemap import Basemap
from .io import ScaleIO
from .calc import *
#from scale.proj import *
import gradsio
import sys

#try:
#    from mpi4py import MPI
#except ImportError:
#    pass


__all__ = ['convert']


config = {
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
'windrot': True,
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
    for key, value in list(kwargs.items()):
        if key in config:
            if value == '-':
                pass
            elif key == 'proj':
                for key2, value2 in list(value.items()):
                    if key2 in config[key]:
                        config[key][key2] = value2
                    else:
                        raise KeyError("'{0:s}' is not a key of config['{1:s}'].".format(key2, key))
            else:
                config[key] = value
        else:
            raise KeyError("'{0:s}' is not a configuration key.".format(key))


def convert_sub(sio, bmap, topo, ftype, lprate, varout_3d, varout_2d, necessary, it, tskip_a, dryrun=False, myrankmsg=''):
    """
    """
    X = {}
    t0_ext = None
    X['topo'] = topo

    if ftype == 'restart':
        for ivar, ivarf in ('rho', 'DENS'), ('rhot', 'RHOT'), \
                           ('qv', 'QV'), ('qc', 'QC'), ('qr', 'QR'), ('qi', 'QI'), ('qs', 'QS'), ('qg', 'QG'):
            if ivar in varout_3d or ivar in necessary:
                X[ivar] = sio.readvar(ivarf, t=it)
        for ivar, ivarf in ('rain', 'SFLX_rain'), ('snow', 'SFLX_snow'), ('glon', 'lon'), ('glat', 'lat'):
            if ivar in varout_2d or ivar in necessary:
                X[ivar] = sio.readvar(ivarf, t=it)

        if 'u' in varout_3d or 'v' in varout_3d:
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

        if 'dbz' in varout_3d or 'max_dbz' in varout_2d:
            if not dryrun:
                print(myrankmsg, 'Calculate: dbz, max_dbz')
            X['dbz'], X['max_dbz'] = calc_ref(sio, t=it, dryrun=dryrun)

        if not dryrun:
            print(myrankmsg, 'Calculate: z')
        X['z'], height_h = calc_height(sio, topo=X['topo'], dryrun=dryrun)

        if 'rhosfc' in varout_2d or 'psfc' in varout_2d:
            if not dryrun:
                print(myrankmsg, 'Calculate: rhosfc, psfc')
                X['rhosfc'], X['psfc'] = calc_rhosfc_psfc(sio, rho=X['rho'], pres=X['p'], height=X['z'], topo=X['topo'], t=it)

        if 'slp' in varout_2d or (config['extrap'] and (config['vcoor'] == 'z' or config['vcoor'] == 'p')):
            if not dryrun:
                print(myrankmsg, 'Calculate smoothed lowest-level surface temperature extrapolated from the free atmosphere')
                t0_ext = extrap_z_t0(sio, X['tk'], lprate=lprate, height=X['z'], t=it)

            if 'slp' in varout_2d:
                if not dryrun:
                    print(myrankmsg, 'Calculate: slp')
                    X['slp'] = calc_slp(sio, p0=X['p'][0], t0_ext=t0_ext, height=X['z'], 
                                        qv=X['qv'], qhydro=X['qhydro'], lprate=lprate, t=it)

    elif ftype == 'restart_sprd':
        for ivar, ivarf in ('u', 'DENS'), ('v', 'MOMX'), ('w', 'MOMY'), ('tk', 'MOMZ'), ('p', 'RHOT'), \
                           ('qv', 'QV'), ('qc', 'QC'), ('qr', 'QR'), ('qi', 'QI'), ('qs', 'QS'), ('qg', 'QG'):
            if ivar in varout_3d:
                X[ivar] = sio.readvar(ivarf, t=it)

        if not dryrun:
            print(myrankmsg, 'Destagger: u, v, w')
            if 'u' in varout_3d:
                X['u'] = calc_destagger(X['u'], axis=2, first_grd=False)
            if 'v' in varout_3d:
                X['v'] = calc_destagger(X['v'], axis=1, first_grd=False)
            if 'w' in varout_3d:
                X['w'] = calc_destagger(X['w'], axis=0, first_grd=False)
                X['w'][0,:,:] = 0.

    elif ftype == 'history':
        for ivar, ivarf in ('rho', 'DENS'), ('momx', 'MOMX'), ('momy', 'MOMY'), ('momz', 'MOMZ'), ('rhot', 'RHOT'), \
                           ('qv', 'QV'), ('qc', 'QC'), ('qr', 'QR'), ('qi', 'QI'), ('qs', 'QS'), ('qg', 'QG'), ('qhydro', 'QHYD'), \
                           ('u', 'U'), ('v', 'V'), ('w', 'W'), ('tk', 'T'), ('p', 'PRES'), ('theta', 'PT'), ('rh', 'RH'):
            if ivar in varout_3d or ivar in necessary:
                X[ivar] = sio.readvar(ivarf, t=it)

        if 'u' in varout_3d or 'v' in varout_3d:
            if not dryrun:
                print(myrankmsg, 'Calculate: rotate u, v')
                X['u'], X['v'] = calc_rotate_winds(sio, bmap, u=X['u'], v=X['v'], t=it)

        for ivar, ivarf in ('u10', 'U10'), ('v10', 'V10'), ('t2', 'T2'), ('q2', 'Q2'), ('olr', 'OLR'), ('slp', 'MSLP'), \
                           ('sst', 'OCEAN_TEMP'), ('tsfc', 'SFC_TEMP'), ('tsfcocean', 'OCEAN_SFC_TEMP'), ('glon', 'lon'), ('glat', 'lat'):
            if ivar in varout_2d or ivar in necessary:
                X[ivar] = sio.readvar(ivarf, t=it)

        if 'u10' in varout_2d or 'v10' in varout_2d:
            if not dryrun:
                print(myrankmsg, 'Calculate: rotate u10, v10')
                X['u10'], X['v10'] = calc_rotate_winds(sio, bmap, u=X['u10'], v=X['v10'], t=it)

        for ivar, ivarf in ('rain', 'RAIN'), ('snow', 'SNOW'):
            if ivar in varout_2d or ivar in necessary:
                iits = max(it-tskip_a+1, 0)
                if dryrun:
                    for iit in range(iits, it+1):
                        sio.readvar(ivarf, t=iit)
                else:
                    X[ivar] = sio.readvar(ivarf, t=iits)
                    for iit in range(iits+1, it+1):
                        X[ivar] += sio.readvar(ivarf, t=iit)
                    X[ivar] /= (it - iits + 1)

        if 'qhydro' in varout_3d and 'qhydro' not in X:
            if not dryrun:
                print(myrankmsg, 'Calculate: qhydro')
            X['qhydro'] = calc_qhydro(sio, t=it, dryrun=dryrun)
        if 'dbz' in varout_3d or 'max_dbz' in varout_2d:
            if not dryrun:
                print(myrankmsg, 'Calculate: dbz, max_dbz')
            X['dbz'], X['max_dbz'] = calc_ref(sio, t=it, dryrun=dryrun)

        if not dryrun:
            print(myrankmsg, 'Calculate: z')
        X['z'], height_h = calc_height(sio, topo=X['topo'], dryrun=dryrun)
        if not dryrun:
            X['z'] = X['z'].astype('f4')

        if 'slp' in varout_2d and 'slp' not in X or (config['extrap'] and (config['vcoor'] == 'z' or config['vcoor'] == 'p')):
            if not dryrun:
                print(myrankmsg, 'Calculate smoothed lowest-level surface temperature extrapolated from the free atmosphere')
            t0_ext = extrap_z_t0(sio, X['tk'], lprate=lprate, height=X['z'], t=it, dryrun=dryrun)

            if 'slp' in varout_2d and 'slp' not in X:
                if not dryrun:
                    print(myrankmsg, 'Calculate: slp')
                X['slp'] = calc_slp(sio, p0=X['p'][0], t0_ext=t0_ext, height=X['z'],
                                    qv=X['qv'], qhydro=X['qhydro'], lprate=lprate, t=it, dryrun=dryrun)

    elif ftype == 'history_z':
        sys.exit('not done yet...')

    else:
        raise ValueError("ftype = '{0:s}' is not supported. ftype: {'restart', 'restart_sprd', 'history', 'history_z'}".format(ftype))

    return X, t0_ext


def convert(basename, topo=None, gradsfile='out.dat', ctlfile='auto', t=dt.datetime(2000, 1, 1), tint=dt.timedelta(hours=6), comm=None, **kwargs):
    """
    """
    rc(**kwargs)

    # Determine if the mpi4py is used
    #------------
    if comm is None:
        nprocs = 1
        myrank = 0
        myrankmsg = ''
    else:
        nprocs = comm.Get_size()
        myrank = comm.Get_rank()
        print('<< My rank / total processes = {:d} / {:d} >>'.format(myrank, nprocs))
        myrankmsg = '<< Rank {:6d} >> '.format(myrank)

    # Initial settings
    #------------
    if config['ftype'] in ('history', 'history_z'):
        bufsize = 0
    else:
        bufsize = 2
    if config['ftype'] == 'restart_sprd':
        varout_3d = [i for i in config['varout_3d'] if i in var_3d_sprd]
        varout_2d = [i for i in config['varout_2d'] if i in var_2d_sprd]
    else:
        varout_3d = [i for i in config['varout_3d'] if i in var_3d]
        varout_2d = [i for i in config['varout_2d'] if i in var_2d]

    # Initialize the ScaleIO object and get dimensions using the master process,
    # then broadcast the dimensions to other processes
    #------------
    if myrank == 0:
        sio = ScaleIO(basename, cache=True, bufsize=bufsize, verbose=1)

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
        bmap = set_bmap(sio, config['proj'])

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

    if nprocs > 1:
        nx = comm.bcast(nx, root=0)
        ny = comm.bcast(ny, root=0)
        nz = comm.bcast(nz, root=0)
        nzout = comm.bcast(nzout, root=0)
        nt = comm.bcast(nt, root=0)
        bmap = comm.bcast(bmap, root=0)

    # Determine the time frames handled by each process
    #------------
    its_a = config['tstart']
    ite_a = config['tend']
    tskip_a = config['tskip']
    if ite_a == -1:
        ite_a = nt
    nto_a = len(range(its_a, ite_a, tskip_a))

    its = its_a + tskip_a * myrank
    ite = ite_a
    tskip = tskip_a * nprocs
    nto = len(range(its, ite, tskip))

    # Generate the CTL file using the master process
    #------------
    if myrank == 0 and ctlfile is not None:
        print('Generate CTL file')

        if bufsize == 0:
            sliceobj = slice(None)
        else:
            sliceobj = slice(bufsize, -bufsize)
        lons = np.min(sio.lon[sliceobj])
        lone = np.max(sio.lon[sliceobj])
        lats = np.min(sio.lat[sliceobj])
        late = np.max(sio.lat[sliceobj])
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
        else:
            ts = sio.t[0]
            if len(sio.t) > 1:
                tint = (sio.t[1] - sio.t[0]) * tskip_a

        tint_min = int(round(tint.total_seconds() / 60))
        if tint_min == 0:
            tint_min = 1

        if config['proj']['type'] == 'LC':
            if 'basepoint_x' in config['proj'] and config['proj']['basepoint_x'] is None:
                iref = 0.5*(nx+1)
            else:
                iref = config['proj']['basepoint_x'] / float(dx) + 0.5
            if 'basepoint_y' in config['proj'] and config['proj']['basepoint_y'] is None:
                jref = 0.5*(ny+1)
            else:
                jref = config['proj']['basepoint_y'] / float(dy) + 0.5
            pdef = 'pdef {isize:6d} {jsize:6d} lcc {latref:12.6f} {lonref:12.6f} {iref:.1f} {jref:.1f} {Struelat:12.6f} {Ntruelat:12.6f} {slon:12.6f} {dx:12.6f} {dy:12.6f}'.format(
                   isize=nx, jsize=ny,
                   latref=config['proj']['basepoint_lat'], lonref=config['proj']['basepoint_lon'],
                   iref=iref, jref=jref,
                   Struelat=config['proj']['LC_lat1'], Ntruelat=config['proj']['LC_lat2'],
                   slon=config['proj']['basepoint_lon'], dx=dx, dy=dy)
        else:
            raise ValueError('[Error] Unsupport map projection.')

        varstr = ''
        for ivar in varout_3d:
            varstr += "{varname:<12}{nz:6d} 99 {dscr:s}\n".format(varname=ivar, nz=nzout, dscr=var_3d_name[ivar])
        for ivar in varout_2d:
            varstr += "{varname:<12}{nz:6d} 99 {dscr:s}\n".format(varname=ivar, nz=0, dscr=var_2d_name[ivar])

        if ctlfile is 'auto':
            ctlfile0 = '{0:s}.ctl'.format(gradsfile.rsplit('.', 1)[0])
        else:
            ctlfile0 = ctlfile

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
        'nto':    nto_a,
        'ts':     ts.strftime('%H:%MZ%d%b%Y'),
        'tint':   tint_min,
        'pdef':   pdef,
        'nvar':   len(varout_3d) + len(varout_2d),
        'varstr': varstr
        }

        with open(ctlfile0, 'w') as fc:
            fc.write(template.format(**context))

    # Main computation
    #------------
    if nto > 0:
        necessary = []
        if config['vcoor'] == 'z' and config['ftype'] != 'restart_sprd':
            necessary += ['z']
            if config['extrap']:
                necessary += ['p', 'tk', 'theta', 'rho', 'rhot', 'qv', 'qhydro']

        if config['vcoor'] == 'p' and config['ftype'] != 'restart_sprd':
            necessary += ['p']
            if config['extrap']:
                necessary += ['z', 'tk', 'theta', 'rho', 'rhot', 'qv', 'qhydro']
#            if 'qhydro' in varout_3d:
#                necessary += ['qc', 'qr', 'qi', 'qs', 'qg']
            if 'dbz' in varout_3d or 'max_dbz' in varout_2d:
                necessary += ['rho', 'qr', 'qs', 'qg']
        if 'u' in varout_3d or 'v' in varout_3d:
            necessary += ['u', 'v']
        if 'u10' in varout_2d or 'v10' in varout_2d:
            necessary += ['u10', 'v10']
        if config['ftype'] == 'restart':
            if 'psfc' in varout_2d:
                necessary += ['rho', 'p', 'z']

        dummy = np.zeros(1)

        ito = 0
        for it in range(its, ite, tskip):
            it_a = its_a + tskip * ito + tskip_a * myrank
            ito_a = nprocs * ito + myrank

            # read data in 'dryrun' mode, one process follows the previous process sequentially
            ######
            if comm is not None and ito_a >= 1:
                srank = myrank-1
                if srank < 0:
                    srank = nprocs-1
                comm.Recv(dummy, source=srank, tag=10+ito_a-1)
            ######
            if sio is None:
                sio = ScaleIO(basename, cache=True, bufsize=bufsize, verbose=1)
            if config['ftype'] != 'restart_sprd':
                if topo is None:
                    if config['ftype'] == 'restart':
                        topo = sio.readvar('TOPO')
                    elif config['ftype'] == 'history':
                        topo = sio.readvar('topo', t=0)
                elif type(topo) is str :
                    sio_topo = ScaleIO(topo, cache=True, bufsize=2, verbose=1)
                    topo = sio_topo.readvar('TOPO')
                    del sio_topo
            X, t0_ext = convert_sub(sio, bmap, topo, config['ftype'], config['lprate'], varout_3d, varout_2d, necessary, it, tskip_a, dryrun=True)
            ######
            if comm is not None and ito_a+1 < nto_a:
                drank = myrank + 1
                if drank >= nprocs:
                    drank = 0
                comm.Send(dummy, dest=drank, tag=10+ito_a)
            ######

            # do real computation
            X, t0_ext = convert_sub(sio, bmap, topo, config['ftype'], config['lprate'], varout_3d, varout_2d, necessary, it, tskip_a, myrankmsg=myrankmsg)

            Xitp = {}

            # do not consider 'history_z'...
            # require: X['z']
            #  <extrap> X['p'], X['tk'], X['theta'], X['rho'], X['rhot'], X['qv'], X['qhydro']
            if config['vcoor'] == 'z' and config['ftype'] != 'restart_sprd':
                for ivar in X:
#                    if ivar != 'z' and (ivar in varout_3d or ivar in ['p', 'tk', 'theta', 'rho', 'rhot']):
                    if ivar != 'z' and (ivar in varout_3d):
                        print(myrankmsg, 'Vertical interpolation at Z-coordinate: ', ivar)
                        Xitp[ivar] = interp_z(sio, X[ivar], height=X['z'], t=it, extrap=config['extrap'])
                if 'z' in varout_3d:
                    Xitp['z'] = np.empty((len(sio.z), ny, nx), dtype=sio.z.dtype)
                    for ilev in range(len(sio.z)):
                        Xitp['z'][ilev] = sio.z[ilev]

                if config['extrap']:
                    kws = {}
                    kwslist = ''
                    for ivar in ['p', 'tk', 'theta', 'rho', 'rhot']:
                        if ivar in varout_3d:
                            kws[ivar] = Xitp[ivar]
                            kwslist += ivar + ', '
                    print(myrankmsg, ' Calculate extrapolated values under the surface assuming a constant lapse rate: ' + kwslist[0:-2])
                    extrap_z_pt(sio, p0=X['p'][0], t0_ext=t0_ext, height=X['z'],
                                qv=X['qv'], qhydro=X['qhydro'], lprate=config['lprate'], t=it, **kws)


            # do not consider 'history_z'...
            # require: X['p']
            #  <extrap> X['z'], X['tk'], X['theta'], X['rho'], X['rhot'], X['qv'], X['qhydro']
            if config['vcoor'] == 'p' and config['ftype'] != 'restart_sprd':

                for ivar in X:
                    if ivar != 'p' and (ivar in varout_3d):
                        print(myrankmsg, 'Vertical interpolation at P-coordinate: ', ivar)
                        Xitp[ivar] = interp_p(sio, X[ivar], config['plevels'], p=X['p'], t=it, extrap=config['extrap'])

                if 'p' in varout_3d:
                    varshape = list(X[varout_3d[0]].shape)
                    varshape[0] = len(config['plevels'])
                    Xitp['p'] = np.empty(varshape, dtype=X[ivar].dtype)
                    for ilev in range(len(config['plevels'])):
                        Xitp['p'][ilev] = config['plevels'][ilev]

                if config['extrap']:
                    kws = {}
                    kwslist = ''
                    for ivar in ['z', 'tk', 'theta', 'rho', 'rhot']:
                        if ivar in varout_3d:
                            kws[ivar] = Xitp[ivar]
                            kwslist += ivar + ', '
                    print(myrankmsg, 'Calculate extrapolated values under the surface assuming a constant lapse rate: ' + kwslist[0:-2])
                    extrap_p_zt(sio, config['plevels'], p=X['p'], t0_ext=t0_ext, height=X['z'],
                                qv=X['qv'], qhydro=X['qhydro'], lprate=config['lprate'], t=it, **kws)


            sio.freecache()


#            if config['vcoor'] == 'o' or (config['vcoor'] == 'z' and config['ftype'] == 'history') or config['ftype'] == 'restart_sprd': # This assumes Z_INTERP has been done in history files
            if config['vcoor'] == 'o' or config['ftype'] == 'restart_sprd':
                Xout = X
            else:
                Xout = Xitp
            nv3d = 0
            for ivar in varout_3d:
                if ivar not in Xout:
                    raise ValueError("Output variable '" + ivar + "' has not been calculated.")
                nv3d += 1
            X3d = np.empty((nv3d, nzout, ny, nx), dtype='f4')
            iv3d = 0
            for ivar in varout_3d:
                if type(Xout[ivar]) == ma.MaskedArray:
                    X3d[iv3d] = Xout[ivar].filled(fill_value=config['missing'])
                else:
                    X3d[iv3d] = Xout[ivar]
                iv3d += 1

            # always look for 'X' (instead of 'Xout') for 2D variables
            nv2d = 0
            for ivar in varout_2d:
                if ivar not in X:
                    raise ValueError("Output variable '" + ivar + "' has not been calculated.")
                nv2d += 1
            X2d = np.empty((nv2d, ny, nx), dtype='f4')
            iv2d = 0
            for ivar in varout_2d:
                if type(X[ivar]) == ma.MaskedArray:
                    X2d[iv2d] = X[ivar].filled(fill_value=config['missing'])
                else:
                    X2d[iv2d] = X[ivar]
                iv2d += 1

            ######
            if comm is not None and ito_a >= 1:
                srank = myrank-1
                if srank < 0:
                    srank = nprocs-1
                comm.Recv(dummy, source=srank, tag=10+nto_a+ito_a-1)
            ######
            if ito_a == 0:
                f = open(gradsfile, 'wb')
            else:
                f = open(gradsfile, 'r+b')
            for iv in range(nv3d):
                print(myrankmsg, 'Write 3D variable: {:s} [to = {:d}]'.format(varout_3d[iv], ito_a+1))
                gradsio.writegrads(f, X3d[iv], iv+1, nv3d=nv3d, nv2d=nv2d, t=ito_a+1, nx=nx, ny=ny, nz=nzout, nt=nt)
            for iv in range(nv2d):
                print(myrankmsg, 'Write 2D variable: {:s} [to = {:d}]'.format(varout_2d[iv], ito_a+1))
                gradsio.writegrads(f, X2d[iv], nv3d+iv+1, nv3d=nv3d, nv2d=nv2d, t=ito_a+1, nx=nx, ny=ny, nz=nzout, nt=nt)
            f.close()
            ######
            if comm is not None and ito_a+1 < nto_a:
                drank = myrank + 1
                if drank >= nprocs:
                    drank = 0
                comm.Send(dummy, dest=drank, tag=10+nto_a+ito_a)
            ######

            ito += 1

#        if myrank == 0:
#            f.close()

    if sio is not None:
        del sio
