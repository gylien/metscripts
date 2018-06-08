import numpy as np
import datetime as dt
from scale.letkf import letkfout_grads

#--- use mpi4py
from mpi4py import MPI
comm = MPI.COMM_WORLD
#--- do not use mpi4py
#comm = None
#---


vcoor = 'o'
hcoor = 'o'
plevels = [100000., 92500., 85000., 70000., 50000., 30000., 20000., 10000., 5000., 2000., 1000., 500.]
#varout_3d = ['u', 'v', 'w', 'tk', 'qv', 'qc', 'qr', 'qi', 'qs', 'qg', 'qhydro', 'dbz']  # variables for restart file
varout_3d = ['u', 'v', 'w', 'p', 'tk', 'theta', 'z', 'qv', 'qc', 'qr', 'qi', 'qs', 'qg', 'qhydro', 'rh', 'dbz']                                   # variables for history file
#varout_2d = ['topo', 'slp', 'rain', 'snow', 'max_dbz']                                                                     # variables for restart file
varout_2d = ['topo', 'slp', 'u10', 'v10', 't2', 'q2', 'rain', 'snow', 'max_dbz', 'olr', 'tsfc', 'tsfcocean', 'sst', 'glon', 'glat']               # variables for history file
proj = {
'type': 'LC',
'basepoint_lon': 135.523,
'basepoint_lat': 34.823,
'basepoint_x': 60000.0,
'basepoint_y': 60000.0,
'LC_lat1': 32.5,
'LC_lat2': 37.5
}
extrap = True
dlon = 0.5
dlat = 0.5

letkfoutdir = '.'
topofile = 'const/topo/topo'

stime = dt.datetime(2013,  7, 13,  6, 10,  0)
etime = dt.datetime(2013,  7, 13,  6, 10,  0)
tint = dt.timedelta(seconds=30)

outtype = 'fcst'
member = 1

tstart = 0
tend = -1
tskip = 1

sim_read = 8
pnetcdf = False

letkfout_grads(letkfoutdir, topofile=topofile, proj=proj, stime=stime, etime=etime, tint=tint,
               outtype=outtype, member=member,
               vcoor=vcoor, hcoor=hcoor, plevels=plevels, dlon=dlon, dlat=dlat,
               varout_3d=varout_3d, varout_2d=varout_2d, extrap=extrap,
               tstart=tstart, tend=tend, tskip=tskip, comm=comm, sim_read=sim_read, pnetcdf=pnetcdf)
