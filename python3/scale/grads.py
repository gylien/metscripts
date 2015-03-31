import numpy as np
#import numpy.ma as ma
#import datetime as dt
from .io import ScaleIO
from .calc import *
import gradsio


__all__ = ['convert']


missing = -9.99e33
vcoor = 'z'
varlist_3d = ['u', 'v', 'w', 'p', 't', 'theta', 'qv', 'qc', 'qr', 'qi', 'qs', 'qg', 'qhydro']
# rho
# dbz
varlist_2d = ['topo']
#
# max_dbz, slp


def convert(basename, basename_topo=None, ftype='restart', gradsfile='out'):
    """
    """
    sio = ScaleIO(basename)
    nx = sio.dimdef['len_g']['x']
    ny = sio.dimdef['len_g']['y']
    if ftype == 'restart':
        nx -= 4
        ny -= 4
    nz = sio.dimdef['len']['z'][0]
    if sio.t is None:
        nt = 1
    else:
        nt = len(sio.t)
    print('nx =', nx)
    print('ny =', ny)
    print('nz =', nz)
    print('nt =', nt)

    var = {}
    if basename_topo is None:
        var['topo'] = sio.readvar('TOPO')[2:-2,2:-2]
    else:
        sio_topo = ScaleIO(basename_topo)
        var['topo'] = sio_topo.readvar('TOPO')[2:-2,2:-2]
#    var['topo'] = np.tile(var['topo'], (nt, 1, 1))

    f = open(gradsfile + '.dat', 'w+b')

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

            height, height_h = calc_height(sio, topo=var['topo'])

            for ivar in varlist_3d:


                print(ivar)


                var[ivar] = interp_z(sio, var[ivar], height=height)

        for iv, ivar in enumerate(varlist_3d):
            gradsio.writegrads(f, var[ivar].filled(fill_value=missing), iv+1,
                               nv3d=len(varlist_3d), nv2d=len(varlist_2d), t=it+1, nx=nx, ny=ny, nz=nz, nt=nt)
        for iv, ivar in enumerate(varlist_2d):
            gradsio.writegrads(f, var[ivar].filled(fill_value=missing), len(varlist_3d)+iv+1,
                               nv3d=len(varlist_3d), nv2d=len(varlist_2d), t=it+1, nx=nx, ny=ny, nz=nz, nt=nt)

    f.close()










