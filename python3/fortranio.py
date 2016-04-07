import numpy as np
import struct


def fort_seq_read(f, var=None, dtype='f4'):
    if var is not None:
        dtype = var.dtype
    else:
        dtype = np.dtype(dtype)

    buf = f.read(4)
    if not buf: return False
    nbytes = struct.unpack(dtype.byteorder + 'i', buf)[0]

    if var is None:
        res = np.fromfile(f, dtype=dtype, count=nbytes//dtype.itemsize)
    else:
        if nbytes != var.nbytes:
            raise ValueError('Record lengths mismatch. {:d} in the record header; {:d} in the input ndarray. It may be due to the endian mismatch.'.format(nbytes, var.nbytes))
        var[:] = np.fromfile(f, dtype=dtype, count=var.size).reshape(var.shape)
        res = True

    buf = f.read(4)
    if not buf: return False
    nbytes2 = struct.unpack(dtype.byteorder + 'i', buf)[0]
    if nbytes != nbytes2:
        raise ValueError('Record lengths mismatch. {:d} in the record header; {:d} in the record footer.'.format(nbytes, nbytes2))

    return res


def fort_seq_write(f, var):
    f.write(struct.pack(var.dtype.byteorder + 'i', var.nbytes))
    var.tofile(f)
    f.write(struct.pack(var.dtype.byteorder + 'i', var.nbytes))
    return True
