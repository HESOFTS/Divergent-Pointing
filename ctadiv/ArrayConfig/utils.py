import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u

from descartes import PolygonPatch
from shapely.ops import unary_union, polygonize
from shapely.geometry import mapping, LineString, Point


def calc_mean(table, columns):
    """
    Calculate a mean value of columns

    Parameters
    ----------
    table: astropy.table
    columns: str or array
        names of columns

    Returns
    -------
    np.float
    """

    if np.size(columns) == 1:
        columns = [columns]
    
    mean = []
    for column in columns:
        if column in ["p_x", "p_y", "p_z"]:
            mean.append(np.average(table[column], weights=table["d_tel"]))
        else:
            mean.append(np.mean(table[column]))
    return tuple(mean)
    
def deg2rad(table, toDeg=False):
    if toDeg:
        for par in ["az", "alt", "zn"]:
            table[par] = table[par].to(u.deg)
            table[par].info.format = '{:.3f}'

        table["radius"] = convert_radius(table["radius"], table["focal"], toDeg=toDeg)
        table["radius"].info.format = '{:.3f}'

        table["fov"]    = table["fov"].to(u.deg**2)
        table["fov"].info.format = '{:.3f}'
    else:
        for par in ["az", "alt", "zn"]:
            table[par] = table[par].to(u.rad)
            table[par].info.format = '{:.3f}'
        
        table["radius"] = convert_radius(table["radius"], table["focal"], toDeg=toDeg)
        table["radius"].info.format = '{:.3f}'
        
        table["fov"]    = table["fov"].to(u.rad**2)
        table["fov"].info.format = '{:.3f}'
    return table

def convert_radius(radius, focal, toDeg=False):
    if toDeg and radius.unit == u.m:
        temp = np.arctan(np.asarray(radius/focal))
        temp = temp*u.rad
        radius = temp.to(u.deg)
    elif not(toDeg) and radius.unit == u.deg:        
        temp = radius.to(u.rad)
        radius= np.tan(temp.value)*focal
    return radius


def calc_multiplicity(array, plotting=False, ax=None):
    array._table = deg2rad(array.table, toDeg=True)
    
    coord = array.get_pointing_coord(icrs=False)
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

    hfov = []
    for patchsky in result:
         hfov.append(patchsky.area)

    hfov = np.array(hfov)

    # multiplicity associated with each patch
    overlaps = np.array(count_overlaps)
    multiplicity = [[i, hfov[overlaps==i].sum()] for i in set(overlaps)]
    
    if plotting and ax:
        max_multiplicity = array.size_of_array
        cmap = plt.cm.get_cmap('rainbow')
        color_list = cmap(np.linspace(0, 1, max_multiplicity))
        
        for i, pol in enumerate(result):
            colore = overlaps[i]
            vals = np.asarray(mapping(pol)['coordinates'])[0]
            
            ax.add_patch(PolygonPatch(mapping(pol), color=color_list[colore-1]))

        return np.asarray(multiplicity), overlaps, ax
    else:
        return np.asarray(multiplicity)
