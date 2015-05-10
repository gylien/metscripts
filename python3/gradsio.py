import numpy as np
import warnings

def endian_det(endian=None):
    if endian is None:
        return np.dtype('f4')
    elif endian == 'b' or endian == '>':
        return np.dtype('>f4')
    elif endian == 'l' or endian == '<':
        return np.dtype('<f4')
    else:
        raise ValueError("'endian' has to be [None|'b'|'l'|'>'|'<'].")

def readgrads(fo, varid, nv3d=0, nv2d=0, t=1, e=1, nx=1, ny=1, nz=1, nt=1, ne=1, endian=None):
    if varid < 1 or varid > nv3d+nv2d:
        raise ValueError("'varid' is out of range.")
    dtyp = endian_det(endian)

    v_onetime = nx * ny * (nz * nv3d + nv2d)
    vstart = v_onetime * (nt * (e-1) + (t-1))
    if varid <= nv3d:
        vstart += nx * ny * nz * (varid-1)
        vlen = nx * ny * nz
        shape = (nz, ny, nx)
    else:
        vstart += nx * ny * (nz * nv3d + varid-nv3d-1)
        vlen = nx * ny
        shape = (ny, nx)

    fo.seek(4*vstart)
    field = np.fromfile(fo, dtype=dtyp, count=vlen)
    return np.reshape(field, shape)

def writegrads(fo, data, varid, nv3d=0, nv2d=0, t=1, e=1, nx=1, ny=1, nz=1, nt=1, ne=1, endian=None):
    if varid < 1 or varid > nv3d+nv2d:
        raise ValueError("'varid' is out of range.")
    dtyp = endian_det(endian)

    v_onetime = nx * ny * (nz * nv3d + nv2d)
    vstart = v_onetime * (nt * (e-1) + (t-1))
    if varid <= nv3d:
        if data.shape != (nz, ny, nx):
            print(data.shape)
            print((nz,ny,nx))
            raise ValueError("'data' has wrong shape.")
        vstart += nx * ny * nz * (varid-1)
    else:
        if data.shape != (ny, nx):
            raise ValueError("'data' has wrong shape.")
        vstart += nx * ny * (nz * nv3d + varid-nv3d-1)

    if data.dtype != dtyp:
        if data.dtype == np.dtype('>f4') or data.dtype == np.dtype('<f4'):
            data.byteswap(True)
        else:
            warnings.warn('Data type conversion from {0:s} to {1:s}.'.format(str(data.dtype), str(dtyp)))
            data = data.astype(dtyp)

    fo.seek(4*vstart)
    data.tofile(fo)
