import numpy as np
import datetime as dt
import os
#from .io import ScaleIO
from .grads import convert

#try:
#    from mpi4py import MPI
#except ImportError:
#    pass


__all__ = ['letkfout_grads']


def letkfout_grads(letkfoutdir, topofile, proj, stime, etime=None, tint=dt.timedelta(hours=6),
                   outtype='all', member='all', vcoor='o', hcoor='o',
                   plevels=None, dlon=None, dlat=None,
                   varout_3d=None, varout_2d=None,
                   extrap=None, tstart=None, tend=None, tskip=None,
                   comm=None, sim_read=1):
    """
    """
    # Determine if the mpi4py is used
    #------------
    if comm is None:
        nprocs = 1
        myrank = 0
    else:
        nprocs = comm.Get_size()
        myrank = comm.Get_rank()

    if etime is None:
        etime = stime
    if vcoor == 'o':
        vsuffix = 'g'
    elif vcoor == 'z':
        vsuffix = 'gz'
    elif vcoor == 'p':
        vsuffix = 'gp'
    else:
        raise ValueError("'vcoor' must be {'o', 'z', 'p'}.")

    if hcoor == 'all':
        hglist = ['o', 'l']
    elif hcoor in ('o', 'l'):
        hglist = [hcoor]
    elif hasattr(hcoor, '__iter__'):
        hglist = hcoor
    else:
        raise ValueError("'hcoor' must be {'o', 'l', array-like}.")

    hg_suffix = {'o': '', 'l': 'll'}
    hg_gradsargs = {'o': 'gradsfile', 'l': 'gradsfile_ll'}
    hg_ctlargs = {'o': 'ctlfile', 'l': 'ctlfile_ll'}

    mlist = []
    if type(member) == int:
        mlist = ['mean', 'meanf', 'sprd'] + list(range(1, member+1))
    elif member == 'all':
        mlist = ['mean', 'meanf', 'sprd']
    elif member in ('mean', 'meanf', 'sprd'):
        mlist = [member]
    elif hasattr(member, '__iter__'):
        mlist = member
    else:
        raise ValueError("'member' must be {int, 'all', 'mean', 'meanf', 'sprd', array-like}.")

    tlist = []
    if outtype == 'all':
        tlist = ['gues', 'anal', 'fcst']
    elif outtype in ('gues', 'anal', 'fcst'):
        tlist = [outtype]
    elif hasattr(outtype, '__iter__'):
        tlist = outtype
    else:
        raise ValueError("'outtype' must be {'gues', 'anal', 'fcst', array-like}.")

    time = stime
    while time <= etime:
        timef = time.strftime('%Y%m%d%H%M%S')
        timef2 = time.strftime('%Y-%m-%d %H:%M:%S')

        if myrank == 0:
            print('=====================')
            print('[{:s}]'.format(timef2))
            print('=====================')

        for ityp in ['gues', 'anal']:
            if ityp in tlist:
                if member == 'all':
                    im = 1
                    while os.path.isdir("{:s}/{:s}/{:s}/{:04d}".format(letkfoutdir, timef, ityp, im)):
                        mlist.append(im)
                        im += 1

                for im in mlist:
                    if type(im) == int:
                        im = '{:04d}'.format(im)

                    basename = "{:s}/{:s}/{:s}/{:s}/init".format(letkfoutdir, timef, ityp, im)

                    if os.path.isfile(basename + ".pe000000.nc"):
                        kws = {'gradsfile': None, 'ctlfile': None, 'gradsfile_ll': None, 'ctlfile_ll': None}

                        for ihg in hglist:
                            ftype = 'restart'

                            kws[hg_gradsargs[ihg]] = "{:s}/{:s}/{:s}{:s}{:s}/{:s}.grd".format(letkfoutdir, timef, ityp, vsuffix, hg_suffix[ihg], im)
                            if myrank == 0:
                                os.makedirs(os.path.dirname(kws[hg_gradsargs[ihg]]), exist_ok=True)

                            if im == 'mean':
                                kws[hg_ctlargs[ihg]] = "{:s}/ctl/{:s}{:s}{:s}_mean.ctl".format(letkfoutdir, ityp, vsuffix, hg_suffix[ihg])
                            elif im == 'meanf':
                                kws[hg_ctlargs[ihg]] = "{:s}/ctl/{:s}{:s}{:s}_meanf.ctl".format(letkfoutdir, ityp, vsuffix, hg_suffix[ihg])
                            elif im == 'sprd':
                                kws[hg_ctlargs[ihg]] = "{:s}/ctl/{:s}{:s}{:s}_sprd.ctl".format(letkfoutdir, ityp, vsuffix, hg_suffix[ihg])
                                ftype = 'restart_sprd'
                            elif im == '0001':
                                kws[hg_ctlargs[ihg]] = "{:s}/ctl/{:s}{:s}{:s}.ctl".format(letkfoutdir, ityp, vsuffix, hg_suffix[ihg])

                            if kws[hg_ctlargs[ihg]] is not None:
                                if myrank == 0:
                                    os.makedirs(os.path.dirname(kws[hg_ctlargs[ihg]]), exist_ok=True)
                                if os.path.isfile(kws[hg_ctlargs[ihg]]):
                                    kws[hg_ctlargs[ihg]] = None

#                        if myrank == 0:
#                            print(kws)

                        if nprocs > 1:
                            comm.Barrier()
                        if myrank == 0:
                            print()
                            print('* {:s}.pe______.nc'.format(basename))
                            print()
                        convert(basename, topo=topofile, t=time, tint=tint,
                                ftype=ftype, vcoor=vcoor, plevels=plevels, dlon=dlon, dlat=dlat,
                                varout_3d=varout_3d, varout_2d=varout_2d,
                                proj=proj, extrap=extrap, comm=comm, sim_read=sim_read, **kws)

                        if nprocs > 1:
                            comm.Barrier()
                        if myrank == 0:
                            for ihg in hglist:
                                if kws[hg_ctlargs[ihg]] is not None:
                                    with open(kws[hg_ctlargs[ihg]], 'r') as f:
                                        ctltext = f.read()

                                    if im == 'mean' or im == 'meanf' or im == 'sprd':
                                        ctltext = ctltext.replace("{:s}/{:s}{:s}{:s}/{:s}.grd\n".format(timef, ityp, vsuffix, hg_suffix[ihg], im),
                                                                  "%y4%m2%d2%h2%n200/{:s}{:s}{:s}/{:s}.grd\noptions template\n".format(ityp, vsuffix, hg_suffix[ihg], im), 1)
                                        ctltext = ctltext.replace('tdef      1', 'tdef  10000', 1)
                                    elif im == '0001':
                                        edef = ''
                                        ie = 0
                                        for imm in mlist:
                                            if type(imm) == int:
                                                edef += "\n{:04d}".format(imm)
                                                ie += 1
                                        ctltext = ctltext.replace("{:s}/{:s}{:s}{:s}/0001.grd\n".format(timef, ityp, vsuffix, hg_suffix[ihg]),
                                                                  "%y4%m2%d2%h2%n200/{:s}{:s}{:s}/%e.grd\noptions template\n".format(ityp, vsuffix, hg_suffix[ihg]), 1)
                                        ctltext = ctltext.replace('tdef      1', 'tdef  10000', 1)
                                        ctltext = ctltext.replace("\npdef", "\nedef   {:4d} names{:s}\npdef".format(ie, edef), 1)

                                    with open(kws[hg_ctlargs[ihg]], 'w') as f:
                                        f.write(ctltext)

        ityp = 'fcst'
        if ityp in tlist:
            if member == 'all':
                im = 1
                while os.path.isdir("{:s}/{:s}/{:s}/{:04d}".format(letkfoutdir, timef, ityp, im)):
                    mlist.append(im)
                    im += 1

            for im in mlist:
                if type(im) == int:
                    im = '{:04d}'.format(im)

                basename = "{:s}/{:s}/{:s}/{:s}/history".format(letkfoutdir, timef, ityp, im)

                if os.path.isfile(basename + ".pe000000.nc"):
                    ftype = 'history'
                    kws = {'gradsfile': None, 'ctlfile': None, 'gradsfile_ll': None, 'ctlfile_ll': None}

                    for ihg in hglist:
                        kws[hg_gradsargs[ihg]] = "{:s}/{:s}/{:s}{:s}{:s}/{:s}.grd".format(letkfoutdir, timef, ityp, vsuffix, hg_suffix[ihg], im)
                        if myrank == 0:
                            os.makedirs(os.path.dirname(kws[hg_gradsargs[ihg]]), exist_ok=True)

                        if im == 'mean':
                            kws[hg_ctlargs[ihg]] = "{:s}/{:s}/ctl/{:s}{:s}{:s}_mean.ctl".format(letkfoutdir, timef, ityp, vsuffix, hg_suffix[ihg])
                        elif im == '0001':
                            kws[hg_ctlargs[ihg]] = "{:s}/{:s}/ctl/{:s}{:s}{:s}.ctl".format(letkfoutdir, timef, ityp, vsuffix, hg_suffix[ihg])

                        if kws[hg_ctlargs[ihg]] is not None:
                            if myrank == 0:
                                os.makedirs(os.path.dirname(kws[hg_ctlargs[ihg]]), exist_ok=True)
                            if os.path.isfile(kws[hg_ctlargs[ihg]]):
                                kws[hg_ctlargs[ihg]] = None

#                    if myrank == 0:
#                        print(kws)

                    if nprocs > 1:
                        comm.Barrier()
                    if myrank == 0:
                        print()
                        print('* {:s}.pe______.nc'.format(basename))
                        print()
                    convert(basename, t=None,
                            ftype=ftype, vcoor=vcoor, plevels=plevels, dlon=dlon, dlat=dlat,
                            varout_3d=varout_3d, varout_2d=varout_2d,
                            proj=proj, extrap=extrap, tstart=tstart, tend=tend, tskip=tskip, comm=comm, sim_read=sim_read, **kws)

                    if nprocs > 1:
                        comm.Barrier()
                    if myrank == 0:
                        for ihg in hglist:
                            if kws[hg_ctlargs[ihg]] is not None:
                                if im == '0001':
                                    with open(kws[hg_ctlargs[ihg]], 'r') as f:
                                        ctltext = f.read()

                                    edef =''
                                    ie = 0
                                    for imm in mlist:
                                        if type(imm) == int:
                                            edef += "\n{:04d}".format(imm)
                                            ie += 1
                                    ctltext = ctltext.replace("{:s}{:s}{:s}/0001.grd\n".format(ityp, vsuffix, hg_suffix[ihg]),
                                                              "{:s}{:s}{:s}/%e.grd\noptions template\n".format(ityp, vsuffix, hg_suffix[ihg]), 1)
                                    ctltext = ctltext.replace("\npdef", "\nedef   {:4d} names{:s}\npdef".format(ie, edef), 1)

                                    with open(kws[hg_ctlargs[ihg]], 'w') as f:
                                        f.write(ctltext)

        time += tint
