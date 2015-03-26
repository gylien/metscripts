import numpy as np
#import numpy.ma as ma
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.colors as pltcol
import mpl_toolkits.basemap as bm
import scaleio as sio


class ScalePlot:
    timebase = sio.scale_time_0

    # shade/contour
    shade = 'auto'    # None: No shade;   'auto': Yes, auto levels; array-like: Yes, given levels
    shade_cmap = None
    shade_norm = None
    contour = None    # None: No contour; 'auto': Yes, auto levels; array-like: Yes, given levels
    coutour_colors = 'gray'

    # projection/range/map
    proj = 'merc'
    lonmin = None
    lonmax = None
    latmin = None
    latmax = None
    llint = None
    zmin = None
    zmax = None
    zint = None
    pmin = None
    pmax = None
    pint = None
    map_res = 'i'
    map_lwidth = 1.5
    map_lcolor = 'green'

    # label/misc
    title = None
    title_fsize = 16
    cbar_loc = 'bottom'
    cbar_pad = '5%'
    cbar_levels = None
    cbar_unit = None


    def __init__(self, basename):
        self.nproc, self.rootgrps, self.dimdef = sio.scale_open(basename, 'r')
        self.z = sio.scale_read(self.nproc, self.rootgrps, self.dimdef, 'z')[1]
        if self.dimdef['len']['time'][0] is None:
            self.t = None
        else:
            self.t = sio.scale_read(self.nproc, self.rootgrps, self.dimdef, 'time', time='all')[1]
        self.z = sio.scale_read(self.nproc, self.rootgrps, self.dimdef, 'z')[1]
        self.lon = sio.scale_read(self.nproc, self.rootgrps, self.dimdef, 'lon')[1]
        self.lat = sio.scale_read(self.nproc, self.rootgrps, self.dimdef, 'lat')[1]


    def __del__(self):
        sio.scale_close(self.rootgrps)


    def setvar(self, **kwargs):
        for key, value in list(kwargs.items()):
            if hasattr(ScalePlot, key):
                setattr(ScalePlot, key, value)
            else:
                raise KeyError("Class ScalePlot does not support variable '{0:s}'.".format(key))


    def xymap(self, data, ax=None, **kwargs):
        self.setvar(**kwargs)

        lonmin = np.min(self.lon) if ScalePlot.lonmin is None else ScalePlot.lonmin ## may have a bug when the plotting area is across lon=0
        lonmax = np.max(self.lon) if ScalePlot.lonmax is None else ScalePlot.lonmax
        latmin = np.min(self.lat) if ScalePlot.latmin is None else ScalePlot.latmin
        latmax = np.max(self.lat) if ScalePlot.latmax is None else ScalePlot.latmax

        if ScalePlot.llint is None:
            d = max(lonmax - lonmin, latmax - latmin) / 5.
            for llint in (0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1., 2., 5., 10., 20., 30., 45., 60., 90.):
                if llint >= d: break
        else:
            llint = ScalePlot.llint
        lon1 = np.floor(lonmin / llint) * llint
        lon2 = (np.ceil(lonmax / llint) + 1) * llint
        lat1 = np.floor(latmin / llint) * llint
        lat2 = (np.ceil(latmax / llint) + 1) * llint

        bmap = bm.Basemap(ax=ax,
                          projection=ScalePlot.proj,
                          llcrnrlat=latmin, urcrnrlat=latmax,
                          llcrnrlon=lonmin, urcrnrlon=lonmax,
                          resolution=ScalePlot.map_res)
        bmap.drawcoastlines(linewidth=ScalePlot.map_lwidth, color=ScalePlot.map_lcolor)
        bmap.drawmeridians(np.arange(lon1, lon2, llint), labels=[0, 0, 0, 1])
        bmap.drawparallels(np.arange(lat1, lat2, llint), labels=[1, 0, 0, 0])

        m_x, m_y = bmap(self.lon, self.lat)

        csf = None
        cbar = None
        if ScalePlot.shade == 'auto':
            csf = bmap.contourf(m_x, m_y, data, cmap=ScalePlot.shade_cmap, norm=ScalePlot.shade_norm, extend='both')
        elif hasattr(ScalePlot.shade, '__iter__'):  # check if it is array-like
            if ScalePlot.shade_norm is None:
                csf = bmap.contourf(m_x, m_y, data, ScalePlot.shade, cmap=ScalePlot.shade_cmap,
                                    norm=pltcol.BoundaryNorm([-float('inf')] + list(ScalePlot.shade) + [float('inf')],
                                    ScalePlot.shade_cmap.N), extend='both')
            else:
                csf = bmap.contourf(m_x, m_y, data, ScalePlot.shade, cmap=ScalePlot.shade_cmap, norm=ScalePlot.shade_norm, extend='both')
        if csf is not None:
            cbar = bmap.colorbar(csf, location=ScalePlot.cbar_loc, pad=ScalePlot.cbar_pad, format='%g')
            if ScalePlot.cbar_levels is None:
                cbar.set_ticks(csf.levels)
            else:
                cbar.set_ticks(ScalePlot.cbar_levels)
            if ScalePlot.cbar_unit is not None:
                cbar.set_label(ScalePlot.cbar_unit)

        cs = None
        if ScalePlot.contour == 'auto':
            cs = bmap.contour(m_x, m_y, data, colors=ScalePlot.coutour_colors)
        elif hasattr(ScalePlot.contour, '__iter__'):  # check if it is array-like
            cs = bmap.contour(m_x, m_y, data, ScalePlot.contour, colors=ScalePlot.coutour_colors)

        if ScalePlot.title is not None:
            if ax is None:
                plt.title(ScalePlot.title, fontsize=ScalePlot.title_fsize)
            else:
                ax.set_title(ScalePlot.title, fontsize=ScalePlot.title_fsize)

        figobj = {'bmap': bmap, 'csf': csf, 'cs': cs, 'cbar': cbar}
        return figobj


    def zlevel(self, varname, t=1, z=1, dscale=1., ax=None, **kwargs):
        self.setvar(**kwargs)

        if self.t is None:
            it = None
            time = None
        elif type(t) is int:
            if t < 0 or t >= len(self.t):
                raise ValueError("'t' is out of the range.")
            timeobj = sio.scale_gettime(self.t[t])
            it = t
            time = None
        elif type(t) is dt.datetime:
            timeobj = t
            it = None
            time = t - ScalePlot.timebase
        else:
            raise ValueError("type of 't' should be either 'int' or 'datetime.datetime'.")
        if self.t is not None:
            print('time =', timeobj.strftime('%Y-%m-%d %H:%M:%S'))

        dim, var = sio.scale_read(self.nproc, self.rootgrps, self.dimdef, varname, time=time, it=it)

        if 'z' in dim:
            figobj = self.xymap(var[z] * dscale, ax=ax)
        else:
            figobj = self.xymap(var * dscale, ax=ax)
        return figobj
