import numpy as np
import numpy.ma as ma
from fortranio import *
import datetime as dt
import time


endian = '>'
dtype_real = endian + 'f4'
dtype_int = endian + 'i4'

Re = 6370.e3 # Earth radius in (m)


def radarobs_read(filename):
    t0 = time.time()
    data = {}

    f = open(filename, 'rb')

    buf = np.zeros(6, dtype=dtype_real)
    fort_seq_read(f, buf)
    data['time'] = dt.datetime(*buf)

    buf = np.zeros(8, dtype=dtype_real)
    fort_seq_read(f, buf)
    data['radar_lon'] = buf[0]
    data['radar_lat'] = buf[1]
    data['radar_alt'] = buf[2]
    data['beam_wid_h'] = buf[3]
    data['beam_wid_v'] = buf[4]
    data['beam_wid_r'] = buf[5]
    data['lambda'] = buf[6]
    data['undef'] = buf[7]

    buf = np.zeros(4, dtype=dtype_int)
    fort_seq_read(f, buf)
    data['na'] = buf[0]
    data['nr'] = buf[1]
    data['ne'] = buf[2]
    data['nvar'] = buf[3]

    data['azim'] = np.zeros(data['na'], dtype=dtype_real)
    fort_seq_read(f, data['azim'])

    data['radi'] = np.zeros(data['nr'], dtype=dtype_real)
    fort_seq_read(f, data['radi'])

    data['elev'] = np.zeros(data['ne'], dtype=dtype_real)
    fort_seq_read(f, data['elev'])

    buf = np.zeros(1, dtype=dtype_real)
    fort_seq_read(f, buf)
    data['attn_fac'] = buf[0]

    buf = np.zeros((data['ne'], data['nr'], data['na']), dtype=dtype_real)
    for ie in xrange(data['ne']):
        fort_seq_read(f, buf[ie])
    data['ref'] = ma.masked_values(buf, data['undef'])

    for ie in xrange(data['ne']):
        fort_seq_read(f, buf[ie])
    data['wind'] = ma.masked_values(buf, data['undef'])

    for ie in xrange(data['ne']):
        fort_seq_read(f, buf[ie])
    data['attn'] = ma.masked_values(buf, data['undef'])

    for ie in xrange(data['ne']):
        fort_seq_read(f, buf[ie])
    data['qc'] = ma.masked_values(buf, data['undef'])

    f.close()

    print "Radar data '{:s}' was read in {:.3f} seconds".format(filename, time.time() - t0)
    return data


def ll_arc_distance(lon0, lat0, arc_dist, az):
    if arc_dist == 0.:
        lon = lon0
        lat = lat0
    else:
        cdist = np.cos(arc_dist / Re)
        sdist = np.sin(arc_dist / Re)
        sinll1 = np.sin(np.deg2rad(lat0))
        cosll1 = np.cos(np.deg2rad(lat0))
        lat = np.rad2deg(np.arcsin(sinll1 * cdist + cosll1 * sdist * np.cos(np.deg2rad(az))))
        lon = lon0 + np.rad2deg(np.arctan2(sdist * np.sin(np.deg2rad(az)), 
                                           cosll1 * cdist - sinll1 * sdist * np.cos(np.deg2rad(az))))
    return lon, lat


def radar_georeference(data):
    Ns = 1.21
    ke = 4. / 3.

    data['r_earth'] = Re

    data['symhgt'] = np.zeros((data['ne'], data['nr']), dtype=dtype_real)

    for ie in xrange(data['ne']):
        for ir in xrange(data['nr']):
            # Two possibilities the second probably more accurate than the
            # first one, but both assume a standard constant value for the
            # refractivity index.
            data['symhgt'][ie,ir] = data['radar_alt'] \
                               + data['radi'][ir] * np.sin(np.deg2rad(data['elev'][ie])) \
                               + (data['radi'][ir] ** 2) / (2. * Ns * Re)
            data['symhgt'][ie,ir] = data['radar_alt'] \
                               + np.sqrt(data['radi'][ir] ** 2 + (ke * Re) ** 2 +
                                         2. * data['radi'][ir] * ke * Re * np.sin(np.deg2rad(data['elev'][ie]))) \
                               - ke * Re

    data['radi_h'] = np.zeros((data['ne'], data['nr']), dtype=dtype_real)
    data['lon'] = np.zeros((data['ne'], data['nr'], data['na']), dtype=dtype_real)
    data['lat'] = np.zeros((data['ne'], data['nr'], data['na']), dtype=dtype_real)
    data['hgt'] = np.zeros((data['ne'], data['nr'], data['na']), dtype=dtype_real)

    for ie in [9]:
#    for ie in xrange(data['ne']):

        print 'ie =', ie

        # The curvature of the radar beam is not taken into account.
        data['radi_h'][ie] = ke * Re * np.arcsin(data['radi'] * np.cos(np.deg2rad(data['elev'][ie])) / (ke * Re)) # Radar horizontal range

        # distance_on_earth=Re*acos(( -radar.radius'.^2+(Re+radar.altitude)^2+(Re+radar.height(:,kk)).^2 )./(2*(Re+radar.altitude)*(Re+radar.height(:,kk))));
        # radar.X(:,:,kk)=distance_on_earth.*sin(az_matrix*d2r);
        # radar.Y(:,:,kk)=distance_on_earth.*cos(az_matrix*d2r);

        for ir in xrange(data['nr']):
            for ia in xrange(data['na']):
                data['lon'][ie,ir,ia], data['lat'][ie,ir,ia] = \
                    ll_arc_distance(data['radar_lon'], data['radar_lat'], data['radi_h'][ie,ir], data['azim'][ia])

    for ia in xrange(data['na']):
        data['hgt'][:,:,ia] = data['symhgt']

    return True
