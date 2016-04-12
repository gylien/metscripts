#!/usr/bin/env python3

from fortranio import *
import sys

if len(sys.argv)-1 < 2:
    raise ValueError("Usage: fort_seq_to_direct inputfile outputfile [dtype]")

ifile = sys.argv[1]
ofile = sys.argv[2]
if len(sys.argv)-1 >= 3:
    dtype = sys.argv[3]
else:
    dtype = 'f4'

f = open(ifile, 'rb')
fo = open(ofile, 'wb')

while True:
    buf = fort_seq_read(f, dtype=dtype)
    if buf is False:
        break

    buf.tofile(fo)

f.close()
fo.close()

