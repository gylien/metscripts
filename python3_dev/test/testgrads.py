import numpy as np
from gradsio import *

f = open('2008020100.grd', 'rb')
v = readgrads(f, 2, nv3d=6, nv2d=6, nx=192, ny=94, nz=64, endian='b')
f.close()

#v = v + 0.

a = np.array(v, dtype='i')

f = open('test.grd', 'wb')
writegrads(f, a, 1, nv3d=1, nx=192, ny=94, nz=64, endian='>')
f.close()
