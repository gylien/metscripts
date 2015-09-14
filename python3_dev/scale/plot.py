import numpy as np
#import numpy.ma as ma
#import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.colors as pltcol
import mpl_toolkits.basemap as bm
#from .io import ScaleIO


__all__ = ['config', 'rc', 'contour_xy']


config = {
### shade/contour
'shade': 'auto',    # None: No shade;   'auto': Yes, auto levels; array-like: Yes, given levels
'shade_cmap': None,
'shade_norm': None,
'contour': None,    # None: No contour; 'auto': Yes, auto levels; array-like: Yes, given levels
'contour_colors': 'gray',
'contour_lwidth': None,
'contour_label': True,
#'contout_label_colors': None,
#'contout_label_fsize': None,

### projection/range/map
'proj': 'merc',
'lonmin': None,
'lonmax': None,
'latmin': None,
'latmax': None,
'llint': None,
'zmin': None,
'zmax': None,
'zint': None,
'pmin': None,
'pmax': None,
'pint': None,
'map_res': 'i',
'map_lwidth': 1.5,
'map_lcolor': 'green',

### slice of data
'slice': slice(None),

### label/misc
'title': None,
'title_fsize': 16,
'cbar_loc': 'bottom',
'cbar_pad': '5%',
'cbar_levels': None,
'cbar_unit': None
}


def rc(**kwargs):
    for key, value in list(kwargs.items()):
        if key in config:
            config[key] = value
        else:
            raise KeyError("'{0:s}' is not a configuration key.".format(key))


def contour_xy(sio, data, ax=None, **kwargs):
    rc(**kwargs)

    lon = sio.lon[config['slice']]
    lat = sio.lat[config['slice']]
    lonmin = np.min(lon) if config['lonmin'] is None else config['lonmin'] ## may have a bug when the plotting area is across lon=0
    lonmax = np.max(lon) if config['lonmax'] is None else config['lonmax']
    latmin = np.min(lat) if config['latmin'] is None else config['latmin']
    latmax = np.max(lat) if config['latmax'] is None else config['latmax']

    if config['llint'] is None:
        d = max(lonmax - lonmin, latmax - latmin) / 5.
        for llint in (0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1., 2., 5., 10., 20., 30., 45., 60., 90.):
            if llint >= d: break
    else:
        llint = config['llint']
    lon1 = np.floor(lonmin / llint) * llint
    lon2 = (np.ceil(lonmax / llint) + 1) * llint
    lat1 = np.floor(latmin / llint) * llint
    lat2 = (np.ceil(latmax / llint) + 1) * llint

    bmap = bm.Basemap(ax=ax,
                      projection=config['proj'],
                      llcrnrlat=latmin, urcrnrlat=latmax,
                      llcrnrlon=lonmin, urcrnrlon=lonmax,
                      resolution=config['map_res'])
    bmap.drawcoastlines(linewidth=config['map_lwidth'], color=config['map_lcolor'])
    bmap.drawmeridians(np.arange(lon1, lon2, llint), labels=[0, 0, 0, 1])
    bmap.drawparallels(np.arange(lat1, lat2, llint), labels=[1, 0, 0, 0])

    m_x, m_y = bmap(lon, lat)

    data_slice = data[config['slice']]

    csf = None
    cbar = None
    if config['shade'] == 'auto':
        csf = bmap.contourf(m_x, m_y, data_slice, cmap=config['shade_cmap'], norm=config['shade_norm'], extend='both')
    elif hasattr(config['shade'], '__iter__'):  # check if it is array-like
        if config['shade_norm'] is None:
            csf = bmap.contourf(m_x, m_y, data_slice, config['shade'], cmap=config['shade_cmap'],
                                norm=pltcol.BoundaryNorm([-float('inf')] + list(config['shade']) + [float('inf')],
                                config['shade_cmap'].N), extend='both')
        else:
            csf = bmap.contourf(m_x, m_y, data_slice, config['shade'], cmap=config['shade_cmap'], norm=config['shade_norm'], extend='both')
    if csf is not None:
        cbar = bmap.colorbar(csf, location=config['cbar_loc'], pad=config['cbar_pad'], format='%g')
        if config['cbar_levels'] is None:
            cbar.set_ticks(csf.levels)
        else:
            cbar.set_ticks(config['cbar_levels'])
        if config['cbar_unit'] is not None:
            cbar.set_label(config['cbar_unit'])

    cs = None
    if config['contour'] == 'auto':
        cs = bmap.contour(m_x, m_y, data_slice, colors=config['contour_colors'], linewidths=config['contour_lwidth'])
    elif hasattr(config['contour'], '__iter__'):  # check if it is array-like
        cs = bmap.contour(m_x, m_y, data_slice, config['contour'], colors=config['contour_colors'], linewidths=config['contour_lwidth'])
    if cs is not None:
        if config['contour_label'] is not None:
            cs.clabel(inline=True, fmt='%g')

    if config['title'] is not None:
        if ax is None:
            plt.title(config['title'], fontsize=config['title_fsize'])
        else:
            ax.set_title(config['title'], fontsize=config['title_fsize'])

    return {'bmap': bmap, 'csf': csf, 'cs': cs, 'cbar': cbar}
