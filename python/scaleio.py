import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
import os.path

dimlist = ['x', 'xh', 'y', 'yh', 'z', 'zh', 'time', 'time1']
dimsublist = ['x', 'xh', 'y', 'yh']
dimorder = [3, 3, 2, 2, 1, 1, 0, 0]
dimoutorder = ['T', 'Z', 'Y', 'X']

def scale_open(basename, mode='r', dim=None, dim_s=None, dim_e=None):
    ncfile = [basename + '.pe{:06d}.nc'.format(0)]
    if not os.path.isfile(ncfile[0]):
        if mode == 'r':
            raise ValueError, '[Error] File does not exist.'
        elif dim is None:
            raise ValueError, '[Error] File does not exist.'
        else:
#            scale_create_new()
            raise ValueError, '[Error] Scale_create has not been supported yet...'

    rootgrp = [Dataset(ncfile[0], mode)]

    dim = dict.fromkeys(dimlist)
    dim_s = dict.fromkeys(dimsublist)
    dim_e = dict.fromkeys(dimsublist)
    for idim in dimlist:
        if idim in rootgrp[0].dimensions:
            dim[idim] = rootgrp[0].variables[idim][:]
            if idim in dimsublist:
                dim_s[idim] = [0]
                dim_e[idim] = [len(rootgrp[0].dimensions[idim])]
        else:
            dim[idim] = None

    ip = 1
    while True:
        ncfile.append(basename + '.pe{:06d}.nc'.format(ip))
        if not os.path.isfile(ncfile[ip]):
            break

        rootgrp.append(Dataset(ncfile[ip], mode))
        for idim in dimlist:
            if idim in dimsublist:
                lent = len(rootgrp[0].dimensions[idim])
                dimt = rootgrp[ip].variables[idim][:]
                if dimt[0] < dim[idim][0]:
                    for i in xrange(ip):
                        dim_s[idim][i] += lent
                        dim_e[idim][i] += lent
                    dim_s[idim].append(0)
                    dim_e[idim].append(lent)
                    dim[idim] = np.concatenate((dimt, dim[idim]))
                elif dimt[-1] > dim[idim][-1]:
                    dim_s[idim].append(max(dim_e[idim]))
                    dim_e[idim].append(dim_s[idim][-1] + lent)
                    dim[idim] = np.concatenate((dim[idim], dimt))
                else:
                    dim_s[idim].append((np.where(dim[idim] == dimt[0]))[0][0])
                    dim_e[idim].append((np.where(dim[idim] == dimt[-1]))[0][0] + 1)
        ip += 1

    nproc = ip
    return nproc, rootgrp, dim, dim_s, dim_e

def scale_close(rootgrp):
    for irg in rootgrp:
        irg.close()

def scale_read(nproc, rootgrp, dim, dim_s, dim_e, varname, time=None, it=None):
    dimtest = rootgrp[0].variables[varname].dimensions
    dimord = [None] * len(dimoutorder)
    dimname = [None] * len(dimoutorder)
    for idim, idimord in zip(dimlist, dimorder):
        if idim in dimtest:
            if dimord[idimord] is None:
                dimord[idimord] = dimtest.index(idim)
                dimname[idimord] = idim
            else:
                raise ValueError, '[Error] Duplicated dimensions.'

    dimshape = []
    for idim in dimname[1:]:
        if idim is not None:
            dimshape.append(len(dim[idim]))
    if dimname[0] is not None:
        if time is not None:
            found = (np.where(dim[dimname[0]] == time))[0]
            if len(found) == 0:
                raise ValueError, '[Error] Cannot find time = ' + str(time) + ' in the file.'
            else:
                it = found[0]
        elif it is None:
            it = 0




#    print dimord
#    print dimname
#    print dimshape
#    print time, it


    print 0


    if hasattr(rootgrp[0].variables[varname][:], 'fill_value'):
        vardata = ma.zeros(dimshape, dtype=rootgrp[0].variables[varname].dtype,
                                     fill_value=rootgrp[0].variables[varname][:].fill_value)
    else:
        vardata = np.zeros(dimshape, dtype=rootgrp[0].variables[varname].dtype)

    for ip in xrange(nproc):
        vartemp = np.transpose(rootgrp[ip].variables[varname], [i for i in dimord if i is not None])
        if dimname[0] is not None:
            vartemp = vartemp[it]

        if dimname[1] is not None and dimname[2] is not None and dimname[3] is not None:
           vardata[:, dim_s[dimname[2]][ip]:dim_e[dimname[2]][ip],
                      dim_s[dimname[3]][ip]:dim_e[dimname[3]][ip]] = vartemp
        elif dimname[1] is None and dimname[2] is not None and dimname[3] is not None:
           vardata[dim_s[dimname[2]][ip]:dim_e[dimname[2]][ip],
                   dim_s[dimname[3]][ip]:dim_e[dimname[3]][ip]] = vartemp
        else:
            raise ValueError, '[Error] Unsupported data dimensions ' + dimname


        print ip+1


    return vardata

