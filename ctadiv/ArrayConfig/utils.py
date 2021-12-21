import astropy.units as u

from astropy.coordinates import SkyCoord, EarthLocation, AltAz, ICRS
from astropy.time import Time

from astroplan import Observer

class CTA_INFO:
    def __init__(self, loc='Roque de los Muchachos', t_obs='2013-11-01T03:00', observer='Roque', verbose=True):
        self._loc = EarthLocation.of_site(loc)
        self._observer = Observer(location=self.loc, name=observer)
        self._t_obs = Time(t_obs)
        self._altaz = AltAz(location=self.loc, obstime=self._t_obs)
        if verbose:
            print("Location        : ", loc)
            print("Observation time: ", t_obs)

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