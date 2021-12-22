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
    
def deg2rad(table, deg=False):
    if deg:
        table["az"]=table["az"].to(u.deg)
        table["alt"]=table["alt"].to(u.deg)
        table["zn"]=table["zn"].to(u.deg)
        table["fov"]=table["fov"].to(u.deg**2)
    else:
        table["az"]=table["az"].to(u.rad)
        table["alt"]=table["alt"].to(u.rad)
        table["zn"]=table["zn"].to(u.rad)
        table["fov"]=table["fov"].to(u.rad**2)
    return table

def group_table(table, group):
    if group:
        tel_group = table.group_by("focal")
        color = ["tab:blue", "tab:orange", "tab:green"]
    else:
        tel_group = table
        color = ["black"]
    return tel_group, color
