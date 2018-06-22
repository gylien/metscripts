import numpy as np
from radartools import *

radarfiles = "/data15/gylien/PAWR/raw/20170716_letkf_v1/pawr_20170716061000"

sfx_ref = ".10000000.dat"
sfx_vr  = ".20000000.dat"
sfx_qc  = "_pawr_qcf.dat"

data = radarobs_read_toshiba(radarfiles+sfx_ref, fname_vr=radarfiles+sfx_vr, fname_qc=radarfiles+sfx_qc)

