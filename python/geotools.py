import numpy as np

Rearth = 6.3712e6 # meter

def distlonlat(lon1, lat1, lon2, lat2):
    lon1r = np.radians(lon1)
    lat1r = np.radians(lat1)
    lon2r = np.radians(lon2)
    lat2r = np.radians(lat2)
    theta = np.arccos(np.cos(lat1r) * np.cos(lat2r) * np.cos(lon1r-lon2r) +
                      np.sin(lat1r) * np.sin(lat2r))
    return Rearth * theta
