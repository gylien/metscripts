import numpy as np
import struct


def fort_seq_read(f, var):
    buf = f.read(4)
    if not buf: return False
    nbytes = struct.unpack(var.dtype.byteorder + 'i', buf)[0]
    if nbytes != var.nbytes:
        raise ValueError('Record lengths mismatch. {:d} in the record header; {:d} in the input ndarray. It may be due to the endian mismatch.'.format(nbytes, var.nbytes))

    var[:] = np.fromfile(f, dtype=var.dtype, count=var.size).reshape(var.shape)

    buf = f.read(4)
    if not buf: return False
    nbytes2 = struct.unpack(var.dtype.byteorder + 'i', buf)[0]
    if nbytes != nbytes2:
        raise ValueError('Record lengths mismatch. {:d} in the record header; {:d} in the record footer.'.format(nbytes, nbytes2))
    return True


def fort_seq_write(f, var):
    f.write(struct.pack(var.dtype.byteorder + 'i', var.nbytes))
    var.tofile(f)
    f.write(struct.pack(var.dtype.byteorder + 'i', var.nbytes))
    return True
