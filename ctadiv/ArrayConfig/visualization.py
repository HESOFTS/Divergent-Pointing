import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt

from ipywidgets import widgets
from IPython.display import display

import copy

import astropy.units as u
from astropy.coordinates import SkyCoord


from astroplan.plots import plot_sky
from astroplan import FixedTarget

from . import utils

import mpl_toolkits.axisartist.angle_helper as angle_helper
from mpl_toolkits.axisartist import SubplotHost, ParasiteAxesAuxTrans
from mpl_toolkits.axisartist.grid_helper_curvelinear import GridHelperCurveLinear

from matplotlib.projections import PolarAxes
from matplotlib.transforms import Affine2D
from astropy.visualization.wcsaxes import SphericalCircle


from descartes import PolygonPatch
from shapely.ops import unary_union, polygonize
from shapely.geometry import mapping, LineString, Point
from ..const import COLORS

def display_1d(table, proj, ax=None, labels=None, **kwargs):
    
    xb = utils.calc_mean(table, proj[0])
    
    ax = plt.figure().add_subplot(111) if ax is None else ax

    for i, [tels, label] in enumerate(zip(table.groups, labels)):
        c = COLORS(i)
        for val in tels[proj[0]]:
            ax.axvline(val, label=label, color=c, **kwargs)
            label='_nolegend_'

    ax.axvline(xb, color="r", label='barycenter', **kwargs)
    ax.set_xlabel("{} [m]".format(proj[0]))
    ax.set_yticks([0, 1])
    ax.legend(frameon=False)

    return ax

def display_2d(table, proj, ax=None, labels=None, **kwargs):
    
    if ax is None:
        ax = plt.figure().add_subplot(111)
        block=True
    else:
        block=False
            
    scale = 1
    
    xb = utils.calc_mean(table, proj[0])
    yb = utils.calc_mean(table, proj[1])
    xbv = utils.calc_mean(table, "p_"+proj[0])
    ybv = utils.calc_mean(table, "p_"+proj[1])

    
    for i, [tels, label] in enumerate(zip(table.groups, labels)):
        xx = tels[proj[0]]
        yy = tels[proj[1]]
        xv = tels["p_"+proj[0]]
        yv = tels["p_"+proj[1]]
        ids = tels["id"]

        s = ax.scatter(xx, yy, label=label, **kwargs)
        ax.quiver(xx, yy, xv, yv, color=s.get_facecolor())

        for i, x, y in zip(ids, xx, yy):
            ax.annotate(i, (x,y))
    ax.scatter(xb, yb, marker='+', label='barycenter', color="r")
    ax.quiver(xb, yb, xbv, ybv, color="r")
    ax.set_xlabel("{} [m]".format(proj[0]))
    ax.set_ylabel("{} [m]".format(proj[1]))

    ax.grid('on')
    ax.axis('equal')
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ax.set_xlim(xlim[0] - 0.25 * np.abs(xlim[0]), xlim[1] + 0.25 * np.abs(xlim[1]))
    ax.set_ylim(ylim[0] - 0.25 * np.abs(ylim[0]), ylim[1] + 0.25 * np.abs(ylim[1]))
    ax.legend(frameon=False)

    return ax

def display_3d(table, proj, ax=None, labels=None, **kwargs):

    ax = plt.figure().add_subplot(111, projection='3d')

    scale = 1

    max_range = []
    for axis in ["x", "y", "z"]:
        max_range.append(table[axis].max() - table[axis].min())
    max_range = max(max_range)

    for i, [tels, label] in enumerate(zip(table.groups, labels)):
        xx = tels["x"]
        yy = tels["y"]
        zz = tels["z"]
        c = COLORS(i)
        ax.quiver(xx, yy, zz, 
                tels["p_x"], tels["p_y"], tels["p_z"],
                length=max_range,
                label=label,
                color=c,
                )

        Xb = scale * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][0].flatten() + scale * (xx.max() + xx.min())
        Yb = scale * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][1].flatten() + scale * (yy.max() + yy.min())
        Zb = scale * max_range * np.mgrid[-0.01:2:2, -0.01:2:2, -0.01:2:2][2].flatten() + scale * (zz.max() + zz.min())
        
        for xb, yb, zb in zip(Xb, Yb, Zb):
            ax.plot([xb], [yb], [zb], 'w')
        
    xx = utils.calc_mean(table, proj[0])
    yy = utils.calc_mean(table, proj[1])
    zz = utils.calc_mean(table, proj[2])
    xbv = utils.calc_mean(table, "p_"+proj[0])
    ybv = utils.calc_mean(table, "p_"+proj[1])
    zbv = utils.calc_mean(table, "p_"+proj[2])

    ax.quiver(xx, yy, zz, 
            xbv, ybv, zbv,
            color="r",
            length=max_range,
            label='barycenter',
            )
     
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.set_zlabel('z [m]')
    ax.legend(frameon=False)
    
    return ax
    
def display_skymap(radec, frame, ax=None, **kwargs):
  
    ax = plt.figure().add_subplot(111, projection='polar') if ax is None else ax

    point = SkyCoord(ra=radec.ra, dec=radec.dec)

    target = FixedTarget(coord=point, name="abc")

    plot_sky(target, frame.observer, frame.t_obs, ax=ax, style_kwargs=kwargs)

    return ax

def skymap_polar(array, group=False, fig=None):

    if fig == None:
        fig = plt.figure() 

    # PolarAxes.PolarTransform takes radian. However, we want our coordinate
    # system in degree
    tr = Affine2D().scale(np.pi/180., 1.).translate(+np.pi/2.,0) + PolarAxes.PolarTransform()

    # polar projection, which involves cycle, and also has limits in
    # its coordinates, needs a special method to find the extremes
    # (min, max of the coordinate within the view).

    # 20, 20 : number of sampling points along x, y direction
    n = 20
    extreme_finder = angle_helper.ExtremeFinderCycle(10, 10,
                                                     lon_cycle=360,
                                                     lat_cycle=None,
                                                     lon_minmax=None,
                                                     lat_minmax=(-90, 90),
                                                     )

    # Find a grid values appropriate for the coordinate (degree,
    # minute, second).
    grid_locator1 = angle_helper.LocatorDMS(12)

    tick_formatter1 = angle_helper.FormatterDMS()
    # And also uses an appropriate formatter.  Note that,the
    # acceptable Locator and Formatter class is a bit different than
    # that of mpl's, and you cannot directly use mpl's Locator and
    # Formatter here (but may be possible in the future).

    grid_helper = GridHelperCurveLinear(tr,
                                        extreme_finder=extreme_finder,
                                        grid_locator1=grid_locator1,
                                        tick_formatter1=tick_formatter1
                                        )
    
    

    ax1 = SubplotHost(fig, 1, 1, 1, grid_helper=grid_helper)

    # make ticklabels of right and top axis visible.
    ax1.axis["right"].major_ticklabels.set_visible(False)
    ax1.axis["top"].major_ticklabels.set_visible(False)
    ax1.axis["bottom"].major_ticklabels.set_visible(True)

    fig.add_subplot(ax1)

    # A parasite axes with given transform
    ax2 = ParasiteAxesAuxTrans(ax1, tr, "equal")
    # note that ax2.transData == tr + ax1.transData
    # Anything you draw in ax2 will match the ticks and grids of ax1.
    ax1.parasites.append(ax2)

    array.convert_unit(toDeg=True)
    tel_group, labels = array.group_by(group)
    for tel_table, label in zip(tel_group.groups, labels):
        s = ax1.scatter(tel_table["az"], tel_table["alt"], label=label,
            s=20, edgecolor="black", transform=ax2.transData, zorder=10)
        
        for tel in tel_table:
            r = SphericalCircle((tel["az"] * u.deg, tel["alt"] * u.deg), tel["radius"] * u.deg, 
                                color=s.get_facecolor()[0], alpha=0.1, transform=ax2.transData)
            ax1.add_patch(r)
            ax2.annotate(tel_table["id"], (tel_table["az"], tel_table["alt"]), fontsize=12, xytext=(4, 4), 
                color="black", textcoords='offset pixels', zorder=10)
        

    ax1.grid(True)
    ax1.set_xlabel("Azimuth [deg]", fontsize=20)
    ax1.set_ylabel("Altitude [deg]", fontsize=20)
    ax1.legend()
    plt.show()

def interactive_polar(array, overwrite=True, group=False):

    if overwrite:
        new_array = array
    else:
        new_array = copy.deepcopy(array)
    fig = plt.figure()

    def update(div=0, az=0, alt=0):

        new_array.divergent_pointing(div, az_mean = az*u.deg, alt_mean = alt*u.deg)
        new_array.convert_unit(toDeg=True)
        plt.cla()
        new_array.skymap_polar(group=group, fig=fig)
        fig.canvas.draw_idle()

    div_s = widgets.FloatLogSlider(value=new_array.div, base=10, min=-4, max =0, step=0.2, description='Divergence')
    az_s = widgets.FloatSlider(value=new_array.pointing["az"].value, min=0, max=360, step=0.01, description='Azumith [deg]')
    alt_s = widgets.FloatSlider(value=new_array.pointing["alt"].value, min=0, max=90, step=0.01, description='Altitude [deg]')
    
    ui = widgets.HBox([div_s, alt_s, az_s])
    out = widgets.interactive_output(update, {'div': div_s, 'az': az_s, 'alt': alt_s})
    display(ui, out)

    return new_array


def multiplicity_plot(array, fig=None):

    if fig == None:
        fig = plt.figure(figsize=(10, 4)) 

    coord = array.pointing_coord(icrs=False)
    polygons = [Point(az, alt).buffer(r) for az, alt, r in zip(coord.az.degree, coord.alt.degree, array.table["radius"])]

    rings = [LineString(list(pol.exterior.coords)) for pol in polygons]
    union = unary_union(rings)
    result = [geom for geom in polygonize(union)]

    count_overlaps = []
    for res in result:
        count_overlaps.append(0)
        for pol in polygons:
            if np.isclose(res.difference(pol).area, 0):
                count_overlaps[-1] +=1
    max_multiplicity = array.size_of_array

    cmap = plt.cm.get_cmap('rainbow')
    color_list = cmap(np.linspace(0, 1, max_multiplicity))
    bounds = np.arange(max_multiplicity + 1) + 1

    hfov = []
    for patchsky in result:
         hfov.append(patchsky.area)

    hfov = np.array(hfov)

    # multiplicity associated with each patch
    overlaps = np.array(count_overlaps)
    average_overlap = np.average(overlaps, weights=hfov)
    variance = np.average((overlaps-average_overlap)**2, weights=hfov)

    gs  = mpl.gridspec.GridSpec(1, 2)

    ax = plt.subplot(gs[0])
    ax_cb = fig.add_axes([0.46,0.1,0.01,0.8])
    ax_mul = plt.subplot(gs[1])

    plt.subplots_adjust(wspace=0.4)

    for i, pol in enumerate(result):
        colore = count_overlaps[i]
        vals = np.asarray(mapping(pol)['coordinates'])[0]
        
        ax.add_patch(PolygonPatch(mapping(pol), color=color_list[colore-1]))

    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    cb1 = mpl.colorbar.ColorbarBase(ax_cb,
                                     norm=norm,
                                     cmap=cmap,
                                     boundaries = bounds,
                                     orientation='vertical',
                                     label='Multiplicity')
    cb1.set_ticks(np.arange(max_multiplicity + 1) + 0.5)
    cb1.set_ticklabels(np.arange(max_multiplicity + 1) + 1)

    ax.set_xlabel("Azimuth [deg]")
    ax.set_ylabel("Altitude [deg]")
    ax.set_xlim(array.pointing["az"].value-20, array.pointing["az"].value+20)
    ax.set_ylim(array.pointing["alt"].value-10, array.pointing["alt"].value+10)
    ax.grid(ls=":", alpha=0.5)
    ax.text(0.9, 0.9, r"Overlap: {:.1f} $\pm$ {:.1f}".format(average_overlap, np.sqrt(variance)), 
            ha="right", transform=ax.transAxes)

    ax_mul.bar(list(set(overlaps)), [hfov[overlaps==i].sum() for i in set(overlaps)])
    ax_mul.text(0.9, 0.9, "Sum = {:.0f}".format(hfov.sum()), ha="right", transform=ax_mul.transAxes)
    ax_mul.set_ylabel('HFOV')
    ax_mul.set_xlabel('Multiplicity')

def interactive_multiplicity(array, overwrite=True):

    if overwrite:
        new_array = array
    else:
        new_array = copy.deepcopy(array)
    
    fig = plt.figure(figsize=(10, 4))

    def update(div=0, az=0, alt=0):

        new_array.divergent_pointing(div, az_mean = az*u.deg, alt_mean = alt*u.deg)
        new_array.convert_unit(toDeg=True)
        plt.cla()
        new_array.multiplicity_plot(fig=fig)
        fig.canvas.draw_idle()

    div_s = widgets.FloatLogSlider(value=new_array.div, base=10, min=-4, max =0, step=0.2, description='Divergence')
    az_s = widgets.FloatSlider(value=new_array.pointing["az"].value, min=0, max=360, step=0.01, description='Azumith [deg]')
    alt_s = widgets.FloatSlider(value=new_array.pointing["alt"].value, min=0, max=90, step=0.01, description='Altitude [deg]')
    
    ui = widgets.HBox([div_s, alt_s, az_s])
    out = widgets.interactive_output(update, {'div': div_s, 'az': az_s, 'alt': alt_s})
    display(ui, out)

    return new_array