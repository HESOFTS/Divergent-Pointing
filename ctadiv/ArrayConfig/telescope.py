import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
import ipywidgets as widgets

from astropy.table import Table
import astropy.table as tab

from . import utils as utils
from . import visualization as visual
from . import pointing

from ..cta import CTA_Info

from astropy.coordinates import SkyCoord

class Telescope:

    _id = 0
    def __init__(self, x, y, z, focal, camera_radius):
        
        Telescope._id += 1
        self.id = Telescope._id
        self.x = x.to(u.m)
        self.y = y.to(u.m)
        self.z = z.to(u.m)
        self.focal = focal.to(u.m)
        self.camera_radius = camera_radius.to(u.m)
        self.alt = u.Quantity(0, u.rad)
        self.az = u.Quantity(0, u.rad)
        self._table = self.make_table()
        
    def point_to_altaz(self, alt, az):
        self.alt = alt.to(u.rad)
        self.az = az.to(u.rad)
        if self.az < 0:
            self.az += 2*np.pi*u.rad
        self.update_table(["alt", "az", "zn", "pointing"])

    @property
    def table(self):
        return self._table

    @property
    def zn(self):
        if self.alt.unit == u.rad:
            return np.pi/2.*u.rad - self.alt
        else:
            return 90*u.deg - self.alt

    @property
    def fov(self):
        """
        Area of the field of view in rad**2
        """
        return np.pi * (self.camera_radius / self.focal)**2*u.rad**2

    @property
    def position(self):
        return np.array([self.x.to(u.m).value, self.y.to(u.m).value, self.z.to(u.m).value]*u.m)

    @property
    def pointing_vector(self):
        # return pointing.alt_az_to_vector(self.alt, self.az)
        return np.array([np.cos(self.alt.to(u.rad))*np.cos(self.az.to(u.rad)),
                         -np.cos(self.alt.to(u.rad))*np.sin(self.az.to(u.rad)),
                         np.sin(self.alt.to(u.rad))])

    def make_table(self):
        properties = [[self.id, self.x.value, self.y.value, self.z.value, 
                       self.az.value, self.alt.value,  self.zn.value, self.focal.value, 
                       self.camera_radius.value, self.fov.value, *tuple(self.pointing_vector),
                       ]]
        
        label = ["id", "x", "y", "z", "az", "alt", "zn", "focal", "radius", "fov", "p_x", "p_y", "p_z"]
        units = ["", u.m, u.m, u.m, u.rad, u.rad, u.rad, u.m, u.m, u.rad**2, "", "", ""]
        dtype = [np.int] + [np.float for i in range(len(units)-1)]
        table = Table(np.asarray(properties, dtype="object"), names=label, units=units, dtype=dtype)
        
        for col in table.columns[4:]:
            table[col].info.format = '{:.3f}'
        
        table.units = "rad"
        return table

    def update_table(self, vals):
        for val in vals:
            if val == "pointing":
                new_val = self.pointing_vector
                for i, v in enumerate(["p_x", "p_y", "p_z"]):
                    self._table[v] = new_val[i]
                    self._table[v].info.format = '{:.3f}'
                
            else:
                self._table[val] = getattr(self, val)
                self._table[val].info.format = '{:.3f}'

    def _convert_unit(self, toDeg=True):
        self._table = utils.deg2rad(self._table, toDeg)
        

class Array:

    def __init__(self, telescope_list, frame=None, **kwargs):
        
        self.telescopes = telescope_list

        self._div = 0
        self._pointing = {"az":0*u.deg, "alt":0*u.deg, "ra": 0*u.deg, "dec": 0*u.deg}

        if frame == None:
            self._frame = CTA_Info(verbose=False, **kwargs)
        else:
            self._frame = frame

        self.__create_table__()

    def __create_table__(self):
        
        table = []
        for tel in self.telescopes:
            table.append(tel.table)
        
        if hasattr(self, "_table"):
            units = self.table.units
        else:
            units = 'rad'

        self._table = tab.vstack(table)

        self._table.add_column(self._dist2tel(), name="d_tel")
        self._table["d_tel"].unit = u.m
        self._table["d_tel"].info.format = "{:.2f}"

        self._table.units = units

    def _convert_unit(self, toDeg=True):
        self._table = utils.deg2rad(self._table, toDeg)
        if toDeg:
            self._table.units = 'deg'
        else:
            self._table.units = 'rad'

    def _dist2tel(self):
        dist = np.zeros(self.size_of_array)
        for i, axis in enumerate(["x", "y", "z"]):
            dist += (self.table[axis] - self.barycenter[i])**2.
        dist = np.sqrt(dist)
        return dist

    @property
    def table(self):
        if hasattr(self._table, "units"):
            if self._table.units == 'deg':
                self._table = utils.deg2rad(self._table, True)
            else:
                self._table = utils.deg2rad(self._table, False)
        return self._table

    @property
    def size_of_array(self):
        return self._table.__len__()

    @property
    def frame(self):
        return self._frame
    
    @property
    def barycenter(self):
        return np.array(utils.calc_mean(self.table, ["x", "y", "z"]))

    @property
    def div(self):
        return self._div

    @property
    def pointing(self):
        return self._pointing

    def update_frame(self, site=None, time=None, delta_t=None, verbose=False):
        self.frame.update(site=site, time=time, delta_t=delta_t, verbose=verbose)
        self.divergent_pointing(self.div, ra=self.pointing["ra"], dec=self.pointing["dec"])

    
    def get_pointing_coord(self, icrs=True):
        return pointing.pointing_coord(self.table, self.frame, icrs=icrs)

    def set_pointing_coord(self, src=None, ra = None, dec = None, alt=None, az = None, unit='deg'):
        
        if type(src) == SkyCoord:
            ra = src.icrs.ra.deg
            dec = src.icrs.dec.deg
            units = 'deg'
            
        if ra is not None and dec is not None:
            src = self.frame.set_source_loc(ra=ra, dec=dec, unit=unit)
            self._pointing["ra"] = src.icrs.ra.value * u.deg
            self._pointing["dec"] = src.icrs.dec.value * u.deg
            self._pointing["alt"] = src.alt.value * u.deg
            self._pointing["az"] = src.az.value * u.deg
        elif alt is not None and az is not None:
            if unit == "deg":
                self._pointing["alt"] = alt * u.deg
                self._pointing["az"] = az * u.deg
            else:
                self._pointing["alt"] = alt * u.rad
                self._pointing["alt"] = az * u.rad
        
        for tel in self.telescopes:
            tel.point_to_altaz(self.pointing["alt"], self.pointing["az"])

        self.__create_table__()
        
    def divergent_pointing(self, div, ra=None, dec = None, alt=None, az=None, unit="deg"):
        """
        Divergent pointing given a parameter div.
        Update pointing of all telescopes of the array

        Parameters
        ----------
        div: float between 0 and 1
        ra(option): float
            source ra 
        dec(option): float
            source dec 
        alt(option): float
            mean alt pointing
        az(option): float
            mean az pointing
        unit(option): string either 'deg' (default) or 'rad'
        """

        self.set_pointing_coord(ra=ra, dec = dec, alt=alt, az=az, unit=unit)

        self._div = div
        
        if div > 1 or div < 0:
            print("[Error] The div value should be between 0 and 1.")
        elif div!=0:
            G = pointing.pointG_position(self.barycenter, self.div, self.pointing["alt"], self.pointing["az"])
            for tel in self.telescopes:
                alt_tel, az_tel = pointing.tel_div_pointing(tel.position, G)
                tel.point_to_altaz(alt_tel*u.rad, az_tel*u.rad)
        
            self.__create_table__()

    def display(self, projection=None, ax=None, group=False, skymap=False, **kwargs):
        """
        Display the array

        Parameters
        ----------
        projection: str
            any combination of 'x', 'y', and 'z'
        ax: `matplotlib.pyplot.axes`
        kwargs: args for `pyplot.scatter` or `pyplot.scatter`

        Returns
        -------
        ax: `matplotlib.pyplot.axes`
        """
        tel_group, labels = self.group_by(group)

        if skymap:
            for i, [table, label] in enumerate(zip(tel_group.groups, labels)):
            
                ax = visual.display_skymap(table, self.frame,  
                                    label=labels[i], ax=ax)
            return ax
        else:
            proj = []
            for axis in ["x", "y", "z"]:
                if axis in projection:
                    proj.append(axis)

        if len(proj) == 1:
            ax = visual.display_1d(tel_group, proj, ax=ax, labels=labels, **kwargs)
        elif len(proj) == 2:
            ax = visual.display_2d(tel_group, proj, ax=ax, labels=labels, **kwargs)
        else:
            ax = visual.display_3d(tel_group, proj, ax=ax, labels=labels, **kwargs)
        
        return ax

    def skymap_polar(self, group=None, fig=None, filename=None):
        return visual.skymap_polar(self, group=group, fig=fig, filename=filename)

    def multiplicity_plot(self, fig=None):
        return visual.multiplicity_plot(self, fig=fig)

    def group_by(self, group = None):
        
        if type(group) == dict:
            groupping = np.zeros(self.size_of_array)
            labels = []
            j = 1
            for key in group.keys():
                labels.append(key)
                for i in group[key]:
                    groupping[i-1] = j
                j+=1
            tel_group = self._table.group_by(np.asarray(groupping))
        elif group:
            tel_group = self._table.group_by("radius")
            labels = ["group_{}".format(i+1) for i in range(len(tel_group.groups))]
        else:
            tel_group = self._table.group_by(np.zeros(self.size_of_array))
            labels = ["_nolegend_"]
        return (tel_group, labels)

    def export_cfg(self, filename=None, verbose=False):
        
        """
        Export cfg file.

        Parameters
        ----------
        filename(option): string
            A default name is 'CTA-ULTRA6-LaPalma-divX-azX-altX.cfg'
        
        verbose(option)

        """

        if filename==None:
            filename = 'CTA-ULTRA6-LaPalma-div{}-az{}-alt{}.cfg'.format(
                str(self.div).replace(".", "_"), 
                str(self.pointing["az"].value).replace(".", "_"), 
                str(self.pointing["alt"].value).replace(".", "_"))
        
        with open(filename, 'w') as f:
            f.write('#ifndef TELESCOPE\n')
            f.write('#  define TELESCOPE 0\n')
            f.write('#endif\n')
            f.write('#if TELESCOPE == 0\n')
            f.write('   TELESCOPE_THETA={:.2f} \n'.format(90 - self.pointing["alt"].value))
            f.write('   TELESCOPE_PHI={:.2f} \n'.format(self.pointing["az"].value))
            f.write('\n% Global and default configuration for things missing in telescope-specific config.\n')
            f.write('#  include <CTA-ULTRA6-LST.cfg>\n')    
            for n, tel in enumerate(self.table, 1):
                zd = 90 - tel['alt']
                f.write('\n#elif TELESCOPE == {:d}\n'.format(n))
                if n <= 4:
                    f.write('#  include <CTA-ULTRA6-LST.cfg>\n')
                else:
                    f.write('#  include <CTA-ULTRA6-MST-NectarCam.cfg>\n')
                f.write('   TELESCOPE_THETA={:.2f}\n'.format(zd))
                f.write('   TELESCOPE_PHI={:.2f}\n'.format(tel["az"]))
            f.write('#else\n')
            f.write('   Error Invalid telescope for CTA-ULTRA6 La Palma configuration.\n')
            f.write('#endif\n')
            f.write('trigger_telescopes = 2 % In contrast to Prod-3 South we apply loose stereo trigger immediately\n')
            f.write('array_trigger = array_trigger_ultra6_diver-test.dat\n')
        
        if verbose:
            with open(filename, 'r') as f:
                for line in f.readlines():
                    print(line)

