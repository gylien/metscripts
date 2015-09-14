import numpy as np
from letkfobsio import *

f = open('1982010100.dat', 'rb')
fw = open('tmp.dat', 'wb')
obs = readobs_all(f, endian='b')
writeobs_all(fw, obs, endian='b')
f.close()
fw.close()

#f = open('1982010100.dat', 'rb')
#fw = open('tmp.dat', 'wb')
#i = 0
#while True:
#    a = readobs(f, endian='b')
#    if not a:
#        break
#    i += 1
#    print i, a
#    writeobs(fw, a, endian='b')
#f.close()
#fw.close()


#f = open('/homes/metogra/gylien/work/miyoshi-read-only/speedy/DATA/obs/RAIN_err_ratio_2p/1982122606.dat', 'rb')
#obs = readobs_all(f, endian='b')
#f.close()

