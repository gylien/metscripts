import numpy as np
import datetime as dt
import os
from .io import ScaleIO
from .grads import convert


__all__ = ['letkfout_grads']


def letkfout_grads(letkfoutdir, topofile, proj, stime, etime=None, tint=dt.timedelta(hours=6),
                   outtype='all', member='all',
                   vcoor='o', plevels='-', varout_3d='-', varout_2d='-', extrap='-', tstart='-', tend='-', tskip='-', threads='-'):
    """
    """
    if etime is None:
        etime = stime
    if vcoor == 'o':
        vsuffix = 'g'
    elif vcoor == 'z':
        vsuffix = 'gz'
    elif vcoor == 'p':
        vsuffix = 'gp'

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

        print('[{:s}]'.format(timef2))

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

                    if os.path.isfile("{:s}/{:s}/{:s}/{:s}/init.pe000000.nc".format(letkfoutdir, timef, ityp, im)):
                        os.makedirs("{:s}/{:s}/{:s}{:s}".format(letkfoutdir, timef, ityp, vsuffix), exist_ok=True)
                        os.makedirs("{:s}/ctl".format(letkfoutdir), exist_ok=True)

                        gradsfile = "{:s}/{:s}/{:s}{:s}/{:s}.grd".format(letkfoutdir, timef, ityp, vsuffix, im)
                        print(gradsfile)

                        ftype = 'restart'
                        if im == 'mean':
                            ctlfile = "{:s}/ctl/{:s}{:s}_mean.ctl".format(letkfoutdir, ityp, vsuffix)
                        elif im == 'meanf':
                            ctlfile = "{:s}/ctl/{:s}{:s}_meanf.ctl".format(letkfoutdir, ityp, vsuffix)
                        elif im == 'sprd':
                            ctlfile = "{:s}/ctl/{:s}{:s}_sprd.ctl".format(letkfoutdir, ityp, vsuffix)
                            ftype = 'restart_sprd'
                        elif im == '0001':
                            ctlfile = "{:s}/ctl/{:s}{:s}.ctl".format(letkfoutdir, ityp, vsuffix)
                        else:
                            ctlfile = None

                        if ctlfile is not None:
                            if os.path.isfile(ctlfile):
                                ctlfile = None  

                        convert("{:s}/{:s}/{:s}/{:s}/init".format(letkfoutdir, timef, ityp, im),
                                topo=topofile, gradsfile=gradsfile, ctlfile=ctlfile, t=time,
                                ftype=ftype, vcoor=vcoor, plevels=plevels, varout_3d=varout_3d, varout_2d=varout_2d, proj=proj, extrap=extrap,
                                threads=threads)

                        if ctlfile is not None:
                            with open(ctlfile, 'r') as f:
                                ctltext = f.read()

                            if im == 'mean' or im == 'meanf' or im == 'sprd':
                                ctltext = ctltext.replace("{:s}/{:s}{:s}/{:s}.grd\n".format(timef, ityp, vsuffix, im),
                                                          "%y4%m2%d2%h2%n200/{:s}{:s}/{:s}.grd\noptions template\n".format(ityp, vsuffix, im), 1)
                                ctltext = ctltext.replace('tdef      1', 'tdef  10000', 1)
                            elif im == '0001':
                                edef = ''
                                ie = 0
                                for imm in mlist:
                                    if type(imm) == int:
                                        edef += "\n{:04d}".format(imm)
                                        ie += 1
                                ctltext = ctltext.replace("{:s}/{:s}{:s}/0001.grd\n".format(timef, ityp, vsuffix),
                                                          "%y4%m2%d2%h2%n200/{:s}{:s}/%e.grd\noptions template\n".format(ityp, vsuffix), 1)
                                ctltext = ctltext.replace('tdef      1', 'tdef  10000', 1)
                                ctltext = ctltext.replace("\npdef", "\nedef   {:4d} names{:s}\npdef".format(ie, edef), 1)

                            with open(ctlfile, 'w') as f:
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

                if os.path.isfile("{:s}/{:s}/{:s}/{:s}/history.pe000000.nc".format(letkfoutdir, timef, ityp, im)):
                    os.makedirs("{:s}/{:s}/{:s}{:s}".format(letkfoutdir, timef, ityp, vsuffix), exist_ok=True)
                    os.makedirs("{:s}/{:s}/ctl".format(letkfoutdir, timef), exist_ok=True)

                    gradsfile = "{:s}/{:s}/{:s}{:s}/{:s}.grd".format(letkfoutdir, timef, ityp, vsuffix, im)
                    print(gradsfile)

                    ftype = 'history'
                    if im == 'mean':
                        ctlfile = "{:s}/{:s}/ctl/{:s}{:s}_mean.ctl".format(letkfoutdir, timef, ityp, vsuffix)
                    elif im == '0001':
                        ctlfile = "{:s}/{:s}/ctl/{:s}{:s}.ctl".format(letkfoutdir, timef, ityp, vsuffix)
                    else:
                        ctlfile = None

                    if ctlfile is not None:
                        if os.path.isfile(ctlfile):
                            ctlfile = None  

                    convert("{:s}/{:s}/{:s}/{:s}/history".format(letkfoutdir, timef, ityp, im),
                    #        topo=topofile, 
                            gradsfile=gradsfile, ctlfile=ctlfile, t=None,
                            ftype=ftype, vcoor=vcoor, plevels=plevels, varout_3d=varout_3d, varout_2d=varout_2d, proj=proj, extrap=extrap,
                            tstart=tstart, tend=tend, tskip=tskip, threads=threads)

                    if ctlfile is not None:
                        if im == '0001':
                            with open(ctlfile, 'r') as f:
                                ctltext = f.read()

                            edef =''
                            ie = 0
                            for imm in mlist:
                                if type(imm) == int:
                                    edef += "\n{:04d}".format(imm)
                                    ie += 1
                            ctltext = ctltext.replace("{:s}{:s}/0001.grd\n".format(ityp, vsuffix),
                                                      "{:s}{:s}/%e.grd\noptions template\n".format(ityp, vsuffix), 1)
                            ctltext = ctltext.replace("\npdef", "\nedef   {:4d} names{:s}\npdef".format(ie, edef), 1)

                            with open(ctlfile, 'w') as f:
                                f.write(ctltext)

        time += tint

