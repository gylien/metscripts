import numpy as np
import datetime as dt
from scale.letkf import letkfout_grads

#--- use mpi4py
from mpi4py import MPI
comm = MPI.COMM_WORLD
#--- do not use mpi4py
#comm = None
#---


sim_read = 12
nprocs = comm.Get_size()
myrank = comm.Get_rank()
if nprocs > sim_read:
    raise ValueError('The maximum number of simultaneous I/O threads is set to ' + str(sim_read) + ', please use nprocs <= ' + str(sim_read))


vcoor = 'o'
hcoor = 'o'
plevels = [100000., 92500., 85000., 70000., 50000., 30000., 20000., 10000., 5000., 2000., 1000., 500.]
varout_3d = ['u', 'v', 'w', 'tk', 'qv', 'qc', 'qr', 'qi', 'qs', 'qg', 'qhydro', 'dbz'] # all available output 3D variables for restart file
varout_2d = ['topo', 'slp', 'rain', 'snow', 'max_dbz']                                 # all available output 2D variables for restart file
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

stime = dt.datetime(2013,  7, 13,  6,  0,  0)
etime = dt.datetime(2013,  7, 13,  6, 10,  0)
tint = dt.timedelta(seconds=30)

outtype = ['anal', 'gues']
member = 1

pnetcdf = False


letkfout_grads(letkfoutdir, topofile=topofile, proj=proj, stime=stime, etime=etime, tint=tint,
               outtype=outtype, member=member,
               vcoor=vcoor, hcoor=hcoor, plevels=plevels, dlon=dlon, dlat=dlat,
               varout_3d=varout_3d, varout_2d=varout_2d, extrap=extrap,
               comm=comm, sim_read=sim_read, pnetcdf=pnetcdf)
