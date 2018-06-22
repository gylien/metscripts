import numpy as np
import numpy.ma as ma
from fortranio import *
import datetime as dt
import time
from c_module.read_toshiba import *


Re = 6370.e3 # Earth radius in (m)


def radarobs_read_toshiba(fname_ref, fname_vr=None, fname_qc=None):
    t0 = time.time()
    data = {}

    hd, az, el, rtdat = read_toshiba(fname_ref)

    try:
        data['time'] = dt.datetime(hd.s_yr, hd.s_mn, hd.s_dy, hd.s_hr, hd.s_mi, hd.s_sc)
    except ValueError:
        data['time'] = None

    data['radar_lon'] = hd.longitude
    data['radar_lat'] = hd.latitude
    data['radar_alt'] = hd.altitude
    data['beam_wid_h'] = hd.beam_wid_h
    data['beam_wid_v'] = hd.beam_wid_v
    data['beam_wid_r'] = float(hd.range_res)
    data['lambda'] = 0.0
    data['undef'] = DMISS

    data['na'] = hd.sector_num
    data['nr'] = hd.range_num
    data['ne'] = hd.total_el_num
    data['nvar'] = 1

    data['azim'] = az[0, 0:data['na']]
    data['radi'] = np.arange(0., float(hd.range_res) * (data['nr']-0.5), float(hd.range_res), dtype='f4')
    data['elev'] = el[0:data['ne'], 0]

    data['attn_fac'] = 0.0

    data['ref'] = ma.masked_values(np.transpose(rtdat[0:data['ne'], 0:data['na'], 0:data['nr']], (0, 2, 1)), data['undef'])

    t1 = time.time()
    print("Radar data '{:s}' was read in {:.3f} seconds".format(fname_ref, t1 - t0))
    t0 = t1

    if fname_vr is not None:
        hd, az, el, rtdat = read_toshiba(fname_vr)
        data['nvar'] += 1
        data['wind'] = ma.masked_values(np.transpose(rtdat[0:data['ne'], 0:data['na'], 0:data['nr']], (0, 2, 1)), data['undef'])

        t1 = time.time()
        print("Radar data '{:s}' was read in {:.3f} seconds".format(fname_vr, t1 - t0))
        t0 = t1

    if fname_qc is not None:
        hd, az, el, rtdat = read_toshiba(fname_qc)
        data['nvar'] += 1
        data['qc'] = ma.masked_values(np.transpose(rtdat[0:data['ne'], 0:data['na'], 0:data['nr']], (0, 2, 1)), data['undef'])

        t1 = time.time()
        print("Radar data '{:s}' was read in {:.3f} seconds".format(fname_qc, t1 - t0))
        t0 = t1

    return data


def radarobs_read(filename, endian=''):
    dtype_real = endian + 'f4'
    dtype_int = endian + 'i4'

    t0 = time.time()
    data = {}

    f = open(filename, 'rb')

    buf = np.zeros(6, dtype=dtype_real)
    fort_seq_read(f, buf)
    try:
        data['time'] = dt.datetime(*buf)
    except ValueError:
        data['time'] = None

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
    for ie in range(data['ne']):
        fort_seq_read(f, buf[ie])
    data['ref'] = ma.masked_values(buf, data['undef'])

    for ie in range(data['ne']):
        fort_seq_read(f, buf[ie])
    data['wind'] = ma.masked_values(buf, data['undef'])

    for ie in range(data['ne']):
        fort_seq_read(f, buf[ie])
    data['qc'] = ma.masked_values(buf, data['undef'])

    for ie in range(data['ne']):
        fort_seq_read(f, buf[ie])
    data['attn'] = ma.masked_values(buf, data['undef'])

    f.close()

    print("Radar data '{:s}' was read in {:.3f} seconds".format(filename, time.time() - t0))
    return data


def dist_ll(lon1, lat1, lon2, lat2):
    return Re * np.arccos(np.sin(np.deg2rad(lat1)) * np.sin(np.deg2rad(lat2)) + 
                          np.cos(np.deg2rad(lat1)) * np.cos(np.deg2rad(lat2)) * np.cos(np.deg2rad(lon2 - lon1)))


def az_ll(lon1, lat1, lon2, lat2):
    return np.rad2deg(np.arctan2(np.sin(np.deg2rad(lon2 - lon1)) * np.cos(np.deg2rad(lat2)),
                                 np.cos(np.deg2rad(lat1)) * np.sin(np.deg2rad(lat2)) - np.sin(np.deg2rad(lat1)) * np.cos(np.deg2rad(lat2)) * np.cos(np.deg2rad(lon2 - lon1))))


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


def radar_georeference(data, lon=None, lat=None, radi_h=None, savedata=False):
    Ns = 1.21
    ke = 4. / 3.

    data['r_earth'] = Re

    data['symhgt'] = np.zeros((data['ne'], data['nr']), dtype='f4')

    for ie in range(data['ne']):
        for ir in range(data['nr']):
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

    data['radi_h'] = np.zeros((data['ne'], data['nr']), dtype='f4')
    data['lon'] = np.zeros((data['ne'], data['nr'], data['na']), dtype='f4')
    data['lat'] = np.zeros((data['ne'], data['nr'], data['na']), dtype='f4')
    data['hgt'] = np.zeros((data['ne'], data['nr'], data['na']), dtype='f4')

    if savedata or (lon is None) or (lat is None) or (radi_h is None):
        for ie in range(data['ne']):

            print('ie =', ie)

            # The curvature of the radar beam is not taken into account.
            data['radi_h'][ie] = ke * Re * np.arcsin(data['radi'] * np.cos(np.deg2rad(data['elev'][ie])) / (ke * Re)) # Radar horizontal range

            for ir in range(data['nr']):
                for ia in range(data['na']):
                    data['lon'][ie,ir,ia], data['lat'][ie,ir,ia] = \
                        ll_arc_distance(data['radar_lon'], data['radar_lat'], data['radi_h'][ie,ir], data['azim'][ia])

        if savedata and (lon is not None) and (lat is not None) and (radi_h is not None):
            np.save(radi_h, data['radi_h'])
            np.save(lon, data['lon'])
            np.save(lat, data['lat'])

    else:
        data['radi_h'] = np.load(radi_h)
        data['lon'] = np.load(lon)
        data['lat'] = np.load(lat)


    for ia in range(data['na']):
        data['hgt'][:,:,ia] = data['symhgt']

    return True
