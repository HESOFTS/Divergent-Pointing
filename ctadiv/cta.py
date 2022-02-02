import numpy as np
import matplotlib.pyplot as plt

from astropy.coordinates import SkyCoord, EarthLocation, AltAz, ICRS
from astropy.time import Time
from astroplan import Observer

import astropy.units as u

from astropy.coordinates import get_moon
from astropy.coordinates import get_sun

    
class CTA_Info:
    def __init__(self, site="North", time = Time.now(), verbose=True):
        self._site = site
        self.site_def(site)
        self._observer = Observer(location=self.loc, name=self.name)
        self._t_obs = Time(time, format='isot', scale='utc')
        self._source = None
        if verbose:
            self.info

    @property
    def info(self):
        print("Observer         : ", self.observer.name)
        print("Location         : ", self.site,  ",", self.loc.to(u.km))
        print("Observation time : ", self.t_obs)

    @property
    def loc(self):
        if self.site.lower() in ('north', 'roque de los muchachos'):
            site_coords = EarthLocation.from_geodetic('342.1184', '28.7606', 2326. * u.meter)
        elif self.site.lower() in ('south', 'paranal'):
            site_coords = EarthLocation.from_geodetic('289.5972', '-24.6253', 2635. * u.meter)
        else:
            raise Warning(f"{site} is not a valid site choice")
        return site_coords
    
    @property
    def t_obs(self):
        return self._t_obs

    @property
    def observer(self):
        return self._observer
    
    @property
    def altaz(self):
        return AltAz(obstime=self.t_obs, location=self.loc)

    @property
    def site(self):
        return self._site

    @property
    def name(self):
        return self._name

    @property
    def delta_t(self):
        return np.linspace(-12, 12, 100)*u.hour

    @property
    def get_sun_loc(self):
        frame = AltAz(obstime=self.t_obs+self.delta_t, location=self.loc)
        return get_sun(self.t_obs+self.delta_t).transform_to(frame)

    @property
    def get_moon_loc(self, time=None):
        frame = AltAz(obstime=self.t_obs+self.delta_t, location=self.loc)
        return get_moon(self.t_obs+self.delta_t).transform_to(frame)

    @property
    def source(self):
        return self._source

    def update(self, site=None, time=None, delta_t=None, verbose=True):
        if site is not None:
            self._site = site
            self.site_def(site)
            self._observer = Observer(location=self.loc, name=self.name)

        elif time is not None:
            self._t_obs = Time(time, format='isot', scale='utc')

        elif delta_t is not None:
            
            if type(delta_t) != u.Quantity:
                print("[Warning] The unit of the input delta_t is assumed to be 'hour'.")
                delta_t = delta_t*u.hour
            
            self._t_obs = Time(self.t_obs+delta_t, format='isot', scale='utc')            
        
        if verbose:
            self.info

    def site_def(self, site):
        if self.site.lower() in ('north', 'roque de los muchachos'):
            self._site = "Roque de los Muchachos"
            self._name = "CTA North"
        elif self.site.lower() in ('south', 'paranal'):
            self._site = "Paranal"
            self._name = "CTA South"
        else:
            raise Warning(f"{site} is not a valid site choice")

    def set_source_loc(self, ra=None, dec=None, timespan=False, unit='deg'):
        if (ra is not None) and (dec is not None):
            source_radec = SkyCoord(ra=ra, dec=dec, frame=ICRS, unit=unit)
        if timespan:
            frame = AltAz(obstime=self.t_obs+self.delta_t, location=self.loc)
        
        self._source = source_radec.transform_to(frame if timespan else self.altaz)
        
        return self._source

    def navigation_plot(self, ra=None, dec=None, unit='deg', **kwargs):
        sun = self.get_sun_loc
        moon = self.get_moon_loc
        plt.plot(self.delta_t, sun.alt, color='r', label='Sun')
        plt.plot(self.delta_t, moon.alt, color=[0.75]*3, ls='--', label='Moon')

        if (ra is not None) and (dec is not None):
            src = self.set_source_loc(ra=ra, dec=dec, timespan=True, unit=unit)
        else:
            src = self.set_source_loc(ra=self.source.icrs.ra, dec=self.source.icrs.dec, timespan=True, unit=unit)

        plt.plot(self.delta_t, src.alt, lw=0.5, alpha=0.5, color="orange")
        plt.scatter(self.delta_t, src.alt,
                    c= src.az, s=8,
                    cmap='viridis',**kwargs)


        plt.fill_between(self.delta_t, 0, 90*u.deg,
                         sun.alt < -0*u.deg, color='0.5', zorder=0)
        plt.fill_between(self.delta_t, 0*u.deg, 90*u.deg,
                         sun.alt < -18*u.deg, color='k', zorder=0)

        plt.colorbar().set_label('Azimuth [deg]')
        plt.legend(loc='upper left')
        plt.xlim(-12*u.hour, 12*u.hour)
        plt.xticks((np.arange(13)*2-12)*u.hour)
        plt.ylim(0*u.deg, 90*u.deg)
        plt.xlabel('Hours from EDT Midnight')
        plt.ylabel('Altitude [deg]')
        plt.show(block=False)
        return plt
