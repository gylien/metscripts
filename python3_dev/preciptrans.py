import numpy as np
import datetime as dt
from gradsio import *
from scipy.stats import norm

ncdf = 200            # number of cdf bins
nTP = 36

ppzero_thres = 0.06   # threshold of no precipitation
mask_thres = 0.35     # threshold of assimilation area wrt. the mask file

# opt_pptrans = 3  # 0: no transformation
#                  # 1: log transformation
#                  # 2: Gaussian transformation with median zero rain
#                  # 3: Gaussian transformation with modified median zero rain

log_trans_tiny = 6.0
gausstail_thres = 0.001

opt_ppobserr = 2 # 0: original obserr form  obs data file
                 # 1: transformed obserr from obs data file
                 # 2: constant obserr
const_ppobserr = 0.5
min_ppobserr = 0.6



def date2tp(datetime):
    if datetime.day <= 10:
        tpm = 1
    elif datetime.day <= 20:
        tpm = 2
    else:
        tpm = 3
    tp = (datetime.month - 1) * 3 + tpm
    return tp


def read_ppcdf(nx, ny, cdfm_dir, cdfo_dir):
    ppcdf_m = np.zeros((nTP,ncdf+1,ny,nx), dtype='f4')
    ppcdf_o = np.zeros((nTP,ncdf+1,ny,nx), dtype='f4')
    ppzero_m = np.zeros((nTP,ny,nx), dtype='f4')
    ppzero_o = np.zeros((nTP,ny,nx), dtype='f4')

    for iTP in range(nTP):
        f = open('{:s}/cdf-{:03d}.dat'.format(cdfm_dir, iTP+1), 'rb')
        ppcdf_m[iTP] = readgrads(f, 1, nv3d=2, nv2d=2, nx=nx, ny=ny, nz=ncdf+1, endian='>')
        ppzero_m[iTP] = readgrads(f, 4, nv3d=2, nv2d=2, nx=nx, ny=ny, nz=ncdf+1, endian='>')
        f.close()
        f = open('{:s}/cdf-{:03d}.dat'.format(cdfo_dir, iTP+1), 'rb')
        ppcdf_o[iTP] = readgrads(f, 1, nv3d=2, nv2d=2, nx=nx, ny=ny, nz=ncdf+1, endian='>')
        ppzero_o[iTP] = readgrads(f, 4, nv3d=2, nv2d=2, nx=nx, ny=ny, nz=ncdf+1, endian='>')
        f.close()

    return (ppcdf_m, ppcdf_o, ppzero_m, ppzero_o)


def compact_tail(pos_cdf):
    if pos_cdf < gausstail_thres:
        return gausstail_thres
    elif pos_cdf > 1.0 - gausstail_thres:
        return 1.0 - gausstail_thres
    else:
        return pos_cdf


#def pptrans(pp, ppcdf, ppzero, method='gauss_cm'):
#    if method not in ['log', 'gauss_cm', 'gauss_bm']:
#        raise ValueError, "'method' has to be ['log'|'gauss_cm'|'gauss_bm']."

#    if ppcdf[0] < -1.0 or ppzero < -1.0:
#        raise ValueError, "Wrong input CDF."

def pptrans_GTcz(pp, ppcdf, ppzero, zerot=ppzero_thres):
    if pp < zerot:
        pos_cdf = ppzero * 0.5
    else:
        if pp < ppcdf[0]:
            pos_cdf = 0.0
        else:
            for b in range(1, ncdf+2):
                if b > ncdf:
                    pos_cdf = 1.0
                    break
                if pp < ppcdf[b]:
                    rr = (pp - ppcdf[b-1]) / (ppcdf[b] - ppcdf[b-1])
                    pos_cdf = ((1.-rr) * (b-1) + rr * b) / ncdf
                    break
    return norm.ppf(compact_tail(pos_cdf))

def pptrans_log(pp, zerot=ppzero_thres, tiny=log_trans_tiny):
    if pp < zerot:
        return np.log(tiny)
    else:
        return np.log(pp + tiny)


