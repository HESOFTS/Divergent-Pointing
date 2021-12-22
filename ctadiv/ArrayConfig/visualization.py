import matplotlib.pyplot as plt
import numpy as np
import mpl_toolkits.axisartist.angle_helper as angle_helper
from mpl_toolkits.axisartist import Subplot
from mpl_toolkits.axisartist import SubplotHost, ParasiteAxesAuxTrans
from mpl_toolkits.axisartist.grid_helper_curvelinear import GridHelperCurveLinear
from matplotlib.projections import PolarAxes
from matplotlib.transforms import Affine2D
import astropy.units as u

from . import utils

def display_1d(table, proj, group=False, ax=None,  **kwargs):
    
    tel_group, color = utils.group_table(table, group)
    xb = utils.calc_mean(table, proj[0])
    
    ax = plt.gca() if ax is None else ax

    for i, [tels, color] in enumerate(zip(tel_group.groups, color)):
        
        if group:
            label = 'group_{}'.format(i)
        else:
            label = '_nolegend_'

        for val in tels[proj[0]]:
            ax.axvline(val, color=color, label=label, **kwargs)
            label='_nolegend_'

    ax.axvline(xb, color="r", label='barycenter', **kwargs)
    ax.set_xlabel("{} [m]".format(proj[0]))
    ax.set_yticks([0, 1])
    ax.legend(frameon=False)
    return ax

def display_2d(table, proj, group=False, ax=None, **kwargs):
    
    tel_group, color = utils.group_table(table, group)

    ax = plt.gca() if ax is None else ax
            
    scale = 1
    
    xb = utils.calc_mean(table, proj[0])
    yb = utils.calc_mean(table, proj[1])
    xbv = utils.calc_mean(table, "p_"+proj[0])
    ybv = utils.calc_mean(table, "p_"+proj[1])

    
    for i, [tels, color] in enumerate(zip(tel_group.groups, color)):
        xx = tels[proj[0]]
        yy = tels[proj[1]]
        xv = tels["p_"+proj[0]]
        yv = tels["p_"+proj[1]]
        
        if group:
            label = 'group_{}'.format(i)
        else:
            label = '_nolegend_'

        ax.scatter(xx, yy, color=color, label=label, **kwargs)
        ax.quiver(xx, yy, xv, yv, color=color)

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

def display_3d(table, proj, group=False, ax=None, **kwargs):

    tel_group, color = utils.group_table(table, group)
    
    ax = plt.figure(figsize=(8, 8)).add_subplot(111, projection='3d')

    scale = 1

    max_range = []
    for axis in ["x", "y", "z"]:
        max_range.append(table[axis].max() - table[axis].min())
    max_range = max(max_range)

    for i, [tels, color] in enumerate(zip(tel_group.groups, color)):
        xx = tels["x"]
        yy = tels["y"]
        zz = tels["z"]

        if group:
            label = 'group_{}'.format(i)
        else:
            label = '_nolegend_'

        ax.quiver(xx, yy, zz, 
                tels["p_x"], tels["p_y"], tels["p_z"],
                color=color,
                length=max_range,
                label=label,
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
    

def display_array_pointing_in_sky(array):
    #TODO: take an array class and plot the pointings FoV in the sky
    # need to handle projection
    pass


def sky_fov(telescope, ax=None):
    """
    Display the telescope FoV in the sky

    Parameters
    ----------
    telescope: `Telescope`
    ax: `matplotlib.pyplot.axes`

    Returns
    -------
    ax: `matplotlib.pyplot.axes`

    """


def polar_stuff(fig, telescope):
    # PolarAxes.PolarTransform takes radian. However, we want our coordinate
    # system in degree
    tr = Affine2D().scale(np.pi/180., 1.).translate(+np.pi/2.,0) + PolarAxes.PolarTransform()

    # polar projection, which involves cycle, and also has limits in
    # its coordinates, needs a special method to find the extremes
    # (min, max of the coordinate within the view).

    # 20, 20 : number of sampling points along x, y direction
    n = 1
    extreme_finder = angle_helper.ExtremeFinderCycle(n, n,
                                                     lon_cycle=360,
                                                     lat_cycle=None,
                                                     lon_minmax=None,
                                                     lat_minmax=(-90, 90),
                                                     )

    grid_locator1 = angle_helper.LocatorDMS(12)
    # Find a grid values appropriate for the coordinate (degree,
    # minute, second).

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
    ax1.axis["right"].major_ticklabels.set_visible(True)
    ax1.axis["top"].major_ticklabels.set_visible(True)

    # let right axis shows ticklabels for 1st coordinate (angle)
    ax1.axis["right"].get_helper().nth_coord_ticks = 0
    # let bottom axis shows ticklabels for 2nd coordinate (radius)
    ax1.axis["bottom"].get_helper().nth_coord_ticks = 1

    fig.add_subplot(ax1)

    # A parasite axes with given transform
    ax2 = ParasiteAxesAuxTrans(ax1, tr, "equal")
    # note that ax2.transData == tr + ax1.transData
    # Anything you draw in ax2 will match the ticks and grids of ax1.
    ax1.parasites.append(ax2)
    # intp = cbook.simple_linear_interpolation
    #ax2.plot(intp(np.array([0, 30]), 50),
    #         intp(np.array([10., 10.]), 50),
    #         linewidth=2.0)



    x = np.rad2deg(telescope.az.value) * np.cos(telescope.alt.value)
    y = np.rad2deg(telescope.alt.value)

    circle = plt.Circle((np.rad2deg(telescope.az.value - np.pi) * np.sin(telescope.alt.value),
                         np.rad2deg(-telescope.alt.value * np.cos((telescope.az.value - np.pi)))),
                        radius=7.7 / 2,
                        color="red",
                        alpha=0.2,
                        )

    circle = plt.Circle((x, y),
                        radius=7.7 / 2,
                        color="red",
                        alpha=0.2,
                        )
    ax1.add_artist(circle)
    # point = ax1.scatter(x, y, c="b", s=20, zorder=10, transform=ax2.transData)
    ax2.annotate(1, (x, y), fontsize=15, xytext=(4, 4), textcoords='offset pixels')

    ax1.set_xlim(-180, 180)
    ax1.set_ylim(0, 90)
    ax1.set_aspect(1.)
    ax1.grid(True, zorder=0)
    ax1.set_xlabel("Azimuth in degrees", fontsize=20)
    ax1.set_ylabel("Zenith in degrees", fontsize=20)

    plt.show()
    return fig



