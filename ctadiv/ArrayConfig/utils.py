import numpy as np

import astropy.units as u

from astropy.coordinates import SkyCoord, EarthLocation, AltAz, ICRS
from astropy.time import Time

from astroplan import Observer

class CTA_INFO:
    def __init__(self, location='Roque de los Muchachos', t_obs='2013-11-01T03:00', observer='Roque', verbose=True):
        self._location = location
        self._loc = EarthLocation.of_site(self.location)

        self._observer = Observer(location=self.loc, name=observer)
        self._t_obs = Time(t_obs)
        self._altaz = AltAz(location=self.loc, obstime=self._t_obs)
        if verbose:
            print("Observer        : ", self.observer.name)
            print("Location        : ", self.location, self.loc, self.loc.to(u.km))
            print("Observation time: ", self.t_obs)

    @property
    def info(self):
        print("Observer        : ", self.observer.name)
        print("Location        : ", self.location, self.loc.to(u.km))
        print("Observation time: ", self.t_obs)

    @property
    def loc(self):
        return self._loc
    
    @property
    def t_obs(self):
        return self._t_obs

    @property
    def observer(self):
        return self._observer
    
    @property
    def altaz(self):
        return self._altaz

    @property
    def location(self):
        return self._location
    


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
    return mean
    
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

def group_table(table, group):
    if group:
        tel_group = table.group_by("radius")
        color = ["tab:blue", "tab:orange", "tab:green"]
    else:
        tel_group = table
        color = ["black"]
    return tel_group, color
