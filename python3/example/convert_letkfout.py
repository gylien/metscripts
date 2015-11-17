import numpy as np
import datetime as dt
from scale.letkf import letkfout_grads
import sys
import os

#--- use mpi4py
from mpi4py import MPI
comm = MPI.COMM_WORLD
#--- do not use mpi4py
#comm = None
#---

vcoor = 'p'
hcoor = 'o'
plevels = [100000., 92500., 85000., 70000., 50000., 30000., 20000., 10000., 5000., 2000., 1000., 500.]

#varout_3d = ['u', 'v', 'w', 'p', 'tk', 'theta', 'rho', 'momx', 'momy', 'momz', 'rhot', 'z', 'qv', 'qc', 'qr', 'qi', 'qs', 'qg', 'qhydro', 'dbz']  # variables for restart file
#varout_3d = ['u', 'v', 'w', 'p', 'tk', 'theta', 'z', 'qv', 'qc', 'qr', 'qi', 'qs', 'qg', 'qhydro', 'rh', 'dbz']                                   # variables for history file
#varout_2d = ['topo', 'rhosfc', 'psfc', 'slp', 'rain', 'snow', 'glon', 'glat']                                                                     # variables for restart file
#varout_2d = ['topo', 'slp', 'u10', 'v10', 't2', 'q2', 'rain', 'snow', 'max_dbz', 'olr', 'tsfc', 'tsfcocean', 'sst', 'glon', 'glat']               # variables for history file
varout_3d = ['u', 'v', 'w', 'p', 'tk', 'z', 'qv', 'qc', 'qr', 'qi', 'qs', 'qg', 'qhydro', 'dbz']
varout_2d = ['topo', 'slp', 'rain', 'snow', 'max_dbz']

proj = {
'type': 'LC',
'basepoint_lon': 140.0,
'basepoint_lat': 30.0,
'basepoint_x': 472000.0,
'basepoint_y': -471000.0,
'LC_lat1': 30.0,
'LC_lat2': 60.0
}
extrap = True
dlon = 0.5
dlat = 0.5

letkfoutdir = '.'
topofile = 'topo/topo'

stime = dt.datetime(2013,  7, 13,  5, 30,  0)
etime = dt.datetime(2013,  7, 13,  6, 00,  0)
tint = dt.timedelta(seconds=1800)

outtype = ['gues', 'anal']
member = 3

tstart = 0
tend = -1
tskip = 1

sim_read = 4

letkfout_grads(letkfoutdir, topofile=topofile, proj=proj, stime=stime, etime=etime, tint=tint,
               outtype=outtype, member=member,
               vcoor=vcoor, hcoor=hcoor, plevels=plevels, dlon=dlon, dlat=dlat,
               varout_3d=varout_3d, varout_2d=varout_2d, extrap=extrap,
               tstart=tstart, tend=tend, tskip=tskip, comm=comm, sim_read=sim_read)
