from ctypes import *
import numpy as np
import numpy.ctypeslib as npct
import os


RDIM = 600
AZDIM = 320
ELDIM = 121
NFNAME = 128
DMISS = -327.68
DNOISE = -327.00

class pawr_header(Structure):
    _fields_ = [("s_yr", c_int), ("s_mn", c_int), ("s_dy", c_int), ("s_hr", c_int), ("s_mi", c_int), ("s_sc", c_int),
                ("e_yr", c_int), ("e_mn", c_int), ("e_dy", c_int), ("e_hr", c_int), ("e_mi", c_int), ("e_sc", c_int),
                ("data_size", c_int),
                ("total_step_num", c_int), ("el_num", c_int), ("total_el_num", c_int),
                ("hit_num", c_int), ("sector_num", c_int), ("range_num", c_int), ("range_res", c_int), ("mesh_size", c_int),
                ("latitude", c_double), ("longitude", c_double), ("altitude", c_double),
                ("start_angle", c_float), ("end_angle", c_float), ("mesh_lsb", c_float), ("mesh_offset", c_float),
                ("tx_freq", c_float), ("tx_power", c_float), ("pulse_len_l", c_float), ("pulse_len_s", c_float),
                ("ant_gain", c_float), ("beam_wid_h", c_float), ("beam_wid_v", c_float),
                ("tx_loss", c_float), ("rx_loss", c_float), ("smin_h", c_float), ("smin_l", c_float),
                ("prf_l", c_float), ("prf_h", c_float), ("zr_b", c_float), ("zr_beta", c_float)]

libdir = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'read_toshiba')
lib_read_toshiba = cdll.LoadLibrary(os.path.join(libdir, "read_toshiba.so"))
lib_read_toshiba.read_toshiba.argtypes = [c_char_p, POINTER(pawr_header), npct.ndpointer(dtype=np.float32), npct.ndpointer(dtype=np.float32), npct.ndpointer(dtype=np.float32)]
lib_read_toshiba.read_toshiba.restype = c_int


def read_toshiba(fname):
    hd = pawr_header()
    az = np.empty((ELDIM, AZDIM), 'f4')
    el = np.empty((ELDIM, AZDIM), 'f4')
    rtdat = np.empty((ELDIM, AZDIM, RDIM), 'f4')

    rc = lib_read_toshiba.read_toshiba(fname.encode(), byref(hd), az, el, rtdat)
    if rc != 0:
        raise IOError('Reading toshiba-format radar data failed; return code = {:d}'.format(rc))

    return hd, az, el, rtdat
