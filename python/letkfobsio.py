import numpy as np
import struct

obsrecords = {
'speedy': {'names':['elm','lon','lat','lev','dat','err'],
           'formats':['i4','f4','f4','f4','f4','f4']},
'gfs':    {'names':['elm','lon','lat','lev','dat','err','typ'],
           'formats':['i4','f4','f4','f4','f4','f4','i4']},
'gfs2':   {'names':['elm','lon','lat','lev','dat','err','typ','tdf','ohx','oqc'],
           'formats':['i4','f4','f4','f4','f4','f4','i4','f4','f4','i4']},
'gfsdep': {'names':['elm','lon','lat','lev','dat','err','typ','tdf','omb','oma'],
           'formats':['i4','f4','f4','f4','f4','f4','i4','f4','f4','f4']}
}

obselems = {
 2819: 'U',
 2820: 'V',
 3073: 'T',
 3074: 'Tv',
 3330: 'Q',
 3331: 'RH',
14593: 'PS',
19999: 'PRC',
99991: 'TCX',
99992: 'TCY',
99993: 'TCP'
}

obstypes = {
 0: 'OTHERS',
 1: 'ADPUPA',
 2: 'AIRCAR',
 3: 'AIRCFT',
 4: 'SATWND',
 5: 'PROFLR',
 6: 'VADWND',
 7: 'SATEMP',
 8: 'ADPSFC',
 9: 'SFCSHP',
10: 'SFCBOG',
11: 'SPSSMI',
12: 'SYNDAT',
13: 'ERS1DA',
14: 'GOESND',
15: 'QKSWND',
16: 'MSONET',
17: 'GPSIPW',
18: 'RASSDA',
19: 'WDSATR',
20: 'ASCATW',
21: 'TMPAPR'
}

def endian_det(endian=None):
    if endian is None:
        return (np.dtype('f4'), 'i')
    elif endian == 'b' or endian == '>':
        return (np.dtype('>f4'), '>i')
    elif endian == 'l' or endian == '<':
        return (np.dtype('<f4'), '<i')
    else:
        raise ValueError, "'endian' has to be [None|'b'|'l'|'>'|'<']."

def readobs(fo, otype='gfs', endian=None):
    try:
        obsrecords[otype]
    except KeyError:
        raise ValueError, "'otype' has to be " + str(obsrecords.keys())
    reclen = len(obsrecords[otype]['names'])
    dtyp, dtypistr = endian_det(endian)

    if not fo.read(4):
        return False
    sglobs = np.fromfile(fo, dtype=dtyp, count=reclen).tolist()
    if not fo.read(4):
        return False

    for i in xrange(reclen):
        if obsrecords[otype]['formats'][i] == 'i4':
            sglobs[i] = int(round(sglobs[i]))
    return sglobs

def readobs_all(fo, otype='gfs', endian=None):
    try:
        obsrecords[otype]
    except KeyError:
        raise ValueError, "'otype' has to be " + str(obsrecords.keys())
    reclen = len(obsrecords[otype]['names'])
    dtyp, dtypistr = endian_det(endian)

    data = []
    fo.seek(0)
    while fo.read(4):
        sglobs = np.fromfile(fo, dtype=dtyp, count=reclen).tolist()
        for i in xrange(reclen):
            if obsrecords[otype]['formats'][i] == 'i4':
                sglobs[i] = int(round(sglobs[i]))
        data.append(tuple(sglobs))
        fo.read(4)

    return np.array(data, dtype=np.dtype(obsrecords[otype]))

def writeobs(fo, sglobs, endian=None):
    dtyp, dtypistr = endian_det(endian)

    sglobs_ndarray = np.array(tuple(sglobs), dtype=dtyp)
    fo.write(struct.pack(dtypistr, sglobs_ndarray.nbytes))
    sglobs_ndarray.tofile(fo)
    fo.write(struct.pack(dtypistr, sglobs_ndarray.nbytes))

def writeobs_all(fo, data, endian=None):
    dtyp, dtypistr = endian_det(endian)

    fo.seek(0)
    for sglobs in data:
        sglobs_ndarray = np.array(tuple(sglobs), dtype=dtyp)
        fo.write(struct.pack(dtypistr, sglobs_ndarray.nbytes))
        sglobs_ndarray.tofile(fo)
        fo.write(struct.pack(dtypistr, sglobs_ndarray.nbytes))
