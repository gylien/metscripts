import numpy as np
#import numpy.ma as ma
import datetime as dt
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as bm
#from cmaptools import *
import scaleio as sio


class scale_plot:
    timebase = sio.scale_time_0


    def __init__(self, basename):
        self.nproc, self.rootgrps, self.dimdef = sio.scale_open(basename, 'r')
        self.z = sio.scale_read(self.nproc, self.rootgrps, self.dimdef, 'z')[1]
        if self.dimdef['len']['time'][0] is None:
            self.t = None
        else:
            self.t = sio.scale_read(self.nproc, self.rootgrps, self.dimdef, 'time', time='all')[1]
        self.lon = sio.scale_read(self.nproc, self.rootgrps, self.dimdef, 'lon')[1]
        self.lat = sio.scale_read(self.nproc, self.rootgrps, self.dimdef, 'lat')[1]


    def __del__(self):
        sio.scale_close(self.rootgrps)


    def xymap(self, data, unit=None, shade='auto', contour=None,
              projection='merc', resolution='i', llint=None, cmap=None, cbar_levels=None,
              title=None, title_fsize=16, display=True, figname=None, dpi=200):
        lonmin = np.min(self.lon) # do not consider area across lon=0
        lonmax = np.max(self.lon)
        latmin = np.min(self.lat)
        latmax = np.max(self.lat)

        if llint is None:
            d = max(lonmax - lonmin, latmax - latmin) / 5.
            for llint in (0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1., 2., 5., 10., 20., 30., 45., 60., 90.):
                if llint >= d: break
        lon1 = np.floor(lonmin / llint) * llint
        lon2 = (np.ceil(lonmax / llint) + 1) * llint
        lat1 = np.floor(latmin / llint) * llint
        lat2 = (np.ceil(latmax / llint) + 1) * llint

        bmap = bm.Basemap(projection=projection,
                          llcrnrlat=latmin, urcrnrlat=latmax,
                          llcrnrlon=lonmin, urcrnrlon=lonmax,
                          resolution=resolution)
        bmap.drawcoastlines(linewidth=1.5, color='green')
#        bmap.drawmapboundary(linewidth=1.0, color='green')
        bmap.drawmeridians(np.arange(lon1, lon2, llint), labels=[0, 0, 0, 1])
        bmap.drawparallels(np.arange(lat1, lat2, llint), labels=[1, 0, 0, 0])

        m_x, m_y = bmap(self.lon, self.lat)

        csf = None
        if shade == 'auto':
            csf = bmap.contourf(m_x, m_y, data, cmap=cmap, extend='both')
        elif hasattr(shade, '__iter__'):  # check if it is array-like
            csf = bmap.contourf(m_x, m_y, data, shade, cmap=cmap, extend='both')
        if csf is not None:
            cbar = bmap.colorbar(csf, location='bottom', pad="5%")
            if cbar_levels is not None:
                cbar.set_ticks(cbar_levels)
            if unit is not None:
                cbar.set_label(unit)

        cs = None
        if contour == 'auto':
            cs = bmap.contour(m_x, m_y, data, colors='gray')
        elif hasattr(contour, '__iter__'):  # check if it is array-like
            cs = bmap.contour(m_x, m_y, data, contour, colors='gray')
#        if cs is not None:
#            plt.clabel(cs, inline=1, fontsize=10)

        if title is not None:
            plt.title(title, fontsize=title_fsize)

        if figname is not None:
            plt.savefig(figname, dpi=dpi)
        if display:
#            plt.show(block=False)
            plt.show()
        plt.clf()

        return True


    def zlevel(self, varname, t=1, z=1, dscale=1., unit=None, shade='auto', contour=None,
               projection='merc', resolution='i', llint=None, cmap=None, cbar_levels=None,
               title=None, title_fsize=16, display=True, figname=None, dpi=200):
        if self.t is None:
            it = None
            time = None
        elif type(t) is int:
            if t < 0 or t >= len(self.t): raise ValueError("'t' is out of the range.")
            timeobj = sio.scale_gettime(self.t[t])
            it = t
            time = None
        elif type(t) is dt.datetime:
            timeobj = t
            it = None
            time = t - scale_plot.timebase
        else:
            raise ValueError("type of 't' should be either 'int' or 'datetime.datetime'.")
        if self.t is not None:
            print 'time =', timeobj.strftime('%Y-%m-%d %H:%M:%S')

        var = sio.scale_read(self.nproc, self.rootgrps, self.dimdef, varname, time=time, it=it)[1]

        self.xymap(var[z] * dscale, unit=unit, shade=shade, contour=contour,
                   projection=projection, resolution=resolution, llint=llint, cmap=cmap, cbar_levels=cbar_levels,
                   title=title, title_fsize=title_fsize, display=display, figname=figname, dpi=dpi)

        return True
