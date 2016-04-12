import numpy as np
import datetime as dt
import os
from .io import ScaleIO
from .grads import convert, conf_default

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
        myrankmsg = ''
    else:
        comm.Barrier()
        nprocs = comm.Get_size()
        myrank = comm.Get_rank()
        myrankmsg = '<< Rank {:6d} >> '.format(myrank)

    # Initial settings
    #------------
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
    elif member in ('mean', 'meanf', 'sprd'):
        mlist = [member]
    elif hasattr(member, '__iter__'):
        mlist = member
    else:
        raise ValueError("'member' must be {int, 'mean', 'meanf', 'sprd', array-like}.")

    tlist = []
    if outtype == 'all':
        tlist = ['gues', 'anal', 'fcst']
    elif outtype in ('gues', 'anal', 'fcst'):
        tlist = [outtype]
    elif hasattr(outtype, '__iter__'):
        tlist = outtype
    else:
        raise ValueError("'outtype' must be {'gues', 'anal', 'fcst', array-like}.")

    # print head messages
    #------------
    if myrank == 0:
        print()
        print('============================================================')
        print('Start time:    {:s}'.format(stime.strftime('%Y-%m-%d %H:%M:%S')))
        print('End time:      {:s}'.format(etime.strftime('%Y-%m-%d %H:%M:%S')))
        print('Time interval: {:.0f} s'.format(tint.total_seconds()))
        print('============================================================')
        print()

    ###### outtype = 'gues', 'anal' ######

    # Determine the total number of files using the master process,
    # then broadcast the dimensions to other processes
    #------------
    n = 0
    basename = []
    kws = []
    time = []
    ftype = []
    ftypeda = []
    im = []

    if myrank == 0:
        ctlfile_done = {
        'gues': dict.fromkeys(('mean', 'meanf', 'sprd', '0001'), False),
        'anal': dict.fromkeys(('mean', 'sprd', '0001'), False)
        }

        itime = stime
        while itime <= etime:
            timef = itime.strftime('%Y%m%d%H%M%S')

            for ityp in ['gues', 'anal']:
                if ityp in tlist:
                    for iim in mlist:
                        if type(iim) == int:
                            iim = '{:04d}'.format(iim)

                        basenametmp = "{:s}/{:s}/{:s}/{:s}/init".format(letkfoutdir, timef, ityp, iim)
                        if os.path.isfile(basenametmp + ".pe000000.nc"):
                            basename.append(basenametmp)
                            kws.append({'gradsfile': None, 'ctlfile': None, 'gradsfile_ll': None, 'ctlfile_ll': None})
                            time.append(itime)
                            if iim == 'sprd':
                                ftype.append('restart_sprd')
                            else:
                                ftype.append('restart')
                            ftypeda.append(ityp)
                            im.append(iim)

                            for ihg in hglist:
                                kws[n][hg_gradsargs[ihg]] = "{:s}/{:s}/{:s}{:s}{:s}/{:s}.grd".format(letkfoutdir, timef, ityp, vsuffix, hg_suffix[ihg], iim)
                                os.makedirs(os.path.dirname(kws[n][hg_gradsargs[ihg]]), exist_ok=True)

                            if iim in ['mean', 'meanf', 'sprd', '0001']:
                                if not ctlfile_done[ityp][iim]:
                                    ctlfile_done[ityp][iim] = True

                                    for ihg in hglist:
                                        if iim == 'mean':
                                            kws[n][hg_ctlargs[ihg]] = "{:s}/ctl/{:s}{:s}{:s}_mean.ctl".format(letkfoutdir, ityp, vsuffix, hg_suffix[ihg])
                                        elif iim == 'meanf':
                                            kws[n][hg_ctlargs[ihg]] = "{:s}/ctl/{:s}{:s}{:s}_meanf.ctl".format(letkfoutdir, ityp, vsuffix, hg_suffix[ihg])
                                        elif iim == 'sprd':
                                            kws[n][hg_ctlargs[ihg]] = "{:s}/ctl/{:s}{:s}{:s}_sprd.ctl".format(letkfoutdir, ityp, vsuffix, hg_suffix[ihg])
                                            ftype[n] = 'restart_sprd'
                                        elif iim == '0001':
                                            kws[n][hg_ctlargs[ihg]] = "{:s}/ctl/{:s}{:s}{:s}.ctl".format(letkfoutdir, ityp, vsuffix, hg_suffix[ihg])

                                        os.makedirs(os.path.dirname(kws[n][hg_ctlargs[ihg]]), exist_ok=True)
                                        if os.path.isfile(kws[n][hg_ctlargs[ihg]]):
                                            kws[n][hg_ctlargs[ihg]] = None

                            n += 1
            itime += tint

    if nprocs > 1:
        n = comm.bcast(n, root=0)
        basename = comm.bcast(basename, root=0)
        kws = comm.bcast(kws, root=0)
        time = comm.bcast(time, root=0)
        ftype = comm.bcast(ftype, root=0)
        ftypeda = comm.bcast(ftypeda, root=0)
        im = comm.bcast(im, root=0)

    # Determine the files handled by each process
    #------------
    nto = 0
    if n > 0:
        if myrank == 0:
            print('------------------------------------------------------------')
            print("Process 'gues' and 'anal' files:")
            for nn in range(n):
                print('{:6d}  {:s}.pe______.nc'.format(nn+1, basename[nn]))
            print('------------------------------------------------------------')
            print()

        nprocs_use = min(nprocs, sim_read)

        its_a = 0
        ite_a = n
        nto_a = len(range(its_a, ite_a))

        if myrank < nprocs_use:
            its = its_a + myrank
            ite = ite_a
            nto = len(range(its, ite, nprocs_use))

    # Main computation
    #------------
    if nprocs > 1:
        comm.Barrier()

    if nto > 0:
        ito = 0
        for it in range(its, ite, nprocs_use):
            it_a = its_a + nprocs_use * ito
            ito_a = nprocs_use * ito + myrank

            print(myrankmsg, '{:s}.pe______.nc'.format(basename[ito_a]))

            convert(basename[ito_a], topo=topofile, t=time[ito_a], tint=tint,
                    ftype=ftype[ito_a], vcoor=vcoor, plevels=plevels, dlon=dlon, dlat=dlat,
                    varout_3d=varout_3d, varout_2d=varout_2d,
                    proj=proj, extrap=extrap, comm=None, sim_read=1, **kws[ito_a])

            for ihg in hglist:
                ifile = kws[ito_a][hg_ctlargs[ihg]]
                if ifile is not None:
                    print(myrankmsg, 'Post-process ctl file: {:s}'.format(ifile))

                    with open(ifile, 'r') as f:
                        ctltext = f.read()

                    if im[ito_a] == 'mean' or im[ito_a] == 'meanf' or im[ito_a] == 'sprd':
                        ctltext = ctltext.replace("{:s}/{:s}{:s}{:s}/{:s}.grd\n".format(time[ito_a].strftime('%Y%m%d%H%M%S'), ftypeda[ito_a], vsuffix, hg_suffix[ihg], im[ito_a]),
                                                  "%y4%m2%d2%h2%n200/{:s}{:s}{:s}/{:s}.grd\noptions template\n".format(ftypeda[ito_a], vsuffix, hg_suffix[ihg], im[ito_a]), 1)
                        ctltext = ctltext.replace('tdef      1', 'tdef  10000', 1)
                    elif im[ito_a] == '0001':
                        edef = ''
                        ie = 0
                        for imm in mlist:
                            if type(imm) == int:
                                edef += "\n{:04d}".format(imm)
                                ie += 1
                        ctltext = ctltext.replace("{:s}/{:s}{:s}{:s}/0001.grd\n".format(time[ito_a].strftime('%Y%m%d%H%M%S'), ftypeda[ito_a], vsuffix, hg_suffix[ihg]),
                                                  "%y4%m2%d2%h2%n200/{:s}{:s}{:s}/%e.grd\noptions template\n".format(ftypeda[ito_a], vsuffix, hg_suffix[ihg]), 1)
                        ctltext = ctltext.replace('tdef      1', 'tdef  10000', 1)
                        ctltext = ctltext.replace("\npdef", "\nedef   {:4d} names{:s}\npdef".format(ie, edef), 1)

                    with open(ifile, 'w') as f:
                        f.write(ctltext)

            ito += 1

    ###### outtype = 'fcst' ######

    # Determine the total number of files using the master process,
    # then broadcast the dimensions to other processes
    #------------
    n = 0
    basename = []
    kws = []
    time = []
    ftype = []
    ftypeda = []
    im = []

    if myrank == 0:
        ctlfile_done = {
        'fcst': dict.fromkeys(('mean', '0001'), False)
        }

        itime = stime
        while itime <= etime:
            timef = itime.strftime('%Y%m%d%H%M%S')

            ityp = 'fcst'
            if ityp in tlist:
                for iim in mlist:
                    if type(iim) == int:
                        iim = '{:04d}'.format(iim)

                    basenametmp = "{:s}/{:s}/{:s}/{:s}/history".format(letkfoutdir, timef, ityp, iim)
                    if os.path.isfile(basenametmp + ".pe000000.nc"):
                        basename.append(basenametmp)
                        kws.append({'gradsfile': None, 'ctlfile': None, 'gradsfile_ll': None, 'ctlfile_ll': None})
                        time.append(itime)
                        ftype.append('history')
                        ftypeda.append(ityp)
                        im.append(iim)

                        for ihg in hglist:
                            kws[n][hg_gradsargs[ihg]] = "{:s}/{:s}/{:s}{:s}{:s}/{:s}.grd".format(letkfoutdir, timef, ityp, vsuffix, hg_suffix[ihg], iim)
                            os.makedirs(os.path.dirname(kws[n][hg_gradsargs[ihg]]), exist_ok=True)

                        if iim in ['mean', 'meanf', 'sprd', '0001']:
                            if not ctlfile_done[ityp][iim]:
                                ctlfile_done[ityp][iim] = True

                                for ihg in hglist:
                                    if iim == 'mean':
                                        kws[n][hg_ctlargs[ihg]] = "{:s}/{:s}/ctl/{:s}{:s}{:s}_mean.ctl".format(letkfoutdir, timef, ityp, vsuffix, hg_suffix[ihg])
                                    elif iim == '0001':
                                        kws[n][hg_ctlargs[ihg]] = "{:s}/{:s}/ctl/{:s}{:s}{:s}.ctl".format(letkfoutdir, timef, ityp, vsuffix, hg_suffix[ihg])

                                    os.makedirs(os.path.dirname(kws[n][hg_ctlargs[ihg]]), exist_ok=True)
                                    if os.path.isfile(kws[n][hg_ctlargs[ihg]]):
                                        kws[n][hg_ctlargs[ihg]] = None

                        n += 1
            itime += tint

    if nprocs > 1:
        n = comm.bcast(n, root=0)
        basename = comm.bcast(basename, root=0)
        kws = comm.bcast(kws, root=0)
        time = comm.bcast(time, root=0)
        ftype = comm.bcast(ftype, root=0)
        ftypeda = comm.bcast(ftypeda, root=0)
        im = comm.bcast(im, root=0)

    # Determine the files handled by each process
    #------------
    nto = 0
    if n > 0:
        if myrank == 0:
            sio = ScaleIO(basename[0], verbose=0)
            if sio.t is None:
                ntf = 1
            else:
                ntf = len(sio.t)
            del sio

            print('------------------------------------------------------------')
            print("Process 'fcst' files:")
            for nn in range(n):
                print('{:6d}  {:s}.pe______.nc'.format(nn+1, basename[nn]))
            print('------------------------------------------------------------')
        else:
            ntf = 0

        if nprocs > 1:
            ntf = comm.bcast(ntf, root=0)

        tstart_r = tstart
        tend_r = tend
        tskip_r = tskip
        if tstart_r is None:
            tstart_r = conf_default['tstart']
        if tend_r is None or tend_r == -1 or tend_r > ntf:
            tend_r = ntf
        if tskip_r is None:
            tskip_r = conf_default['tskip']
        ntfo = len(range(tstart_r, tend_r, tskip_r))

        nsplit = (nprocs-1) // ntfo + 1
        nprocs_L = (nprocs-1) // nsplit + 1
        mycolor = myrank // nprocs_L
        if nprocs > 1:
            commL = comm.Split(mycolor, myrank)
            nprocsL = commL.Get_size()
            myrankL = commL.Get_rank()
        else:
            commL = None
            nprocsL = 1
            myrankL = 0

        sim_read_L = sim_read // min(nsplit, n)

        if myrank == 0:
            print('ntfo =', ntfo)
            print('nsplit =', nsplit)
            print('nprocs_L =', nprocs_L)
            print('sim_read_L =', sim_read_L)
            print('------------------------------------------------------------')
            print()

        its_a = 0
        ite_a = n
        nto_a = len(range(its_a, ite_a))

        if mycolor < nsplit:
            its = its_a + mycolor
            ite = ite_a
            nto = len(range(its, ite, nsplit))

    # Main computation
    #------------
    if nprocs > 1:
        comm.Barrier()

    if nto > 0:
        ito = 0
        for it in range(its, ite, nsplit):
            it_a = its_a + nsplit * ito
            ito_a = nsplit * ito + mycolor

            print(myrankmsg, '{:s}.pe______.nc'.format(basename[ito_a]))

            convert(basename[ito_a], t=None,
                    ftype=ftype[ito_a], vcoor=vcoor, plevels=plevels, dlon=dlon, dlat=dlat,
                    varout_3d=varout_3d, varout_2d=varout_2d,
                    proj=proj, extrap=extrap, tstart=tstart, tend=tend, tskip=tskip,
                    comm=comm, commL=commL, sim_read=sim_read_L, **kws[ito_a])
            if nprocs > 1:
                commL.Barrier()

            if myrankL == 0:
                for ihg in hglist:
                    ifile = kws[ito_a][hg_ctlargs[ihg]]
                    if ifile is not None:
                        print(myrankmsg, 'Post-process ctl file: {:s}'.format(ifile))

                        if im[ito_a] == '0001':
                            with open(ifile, 'r') as f:
                                ctltext = f.read()

                            edef = ''
                            ie = 0
                            for imm in mlist:
                                if type(imm) == int:
                                    edef += "\n{:04d}".format(imm)
                                    ie += 1
                            ctltext = ctltext.replace("{:s}{:s}{:s}/0001.grd\n".format(ftypeda[ito_a], vsuffix, hg_suffix[ihg]),
                                                      "{:s}{:s}{:s}/%e.grd\noptions template\n".format(ftypeda[ito_a], vsuffix, hg_suffix[ihg]), 1)
                            ctltext = ctltext.replace("\npdef", "\nedef   {:4d} names{:s}\npdef".format(ie, edef), 1)

                            with open(ifile, 'w') as f:
                                f.write(ctltext)

            ito += 1
