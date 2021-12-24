import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
import ipywidgets as widgets

from astropy.table import Table
import astropy.table as tab

from . import utils as utils
from . import visualization as visual
from . import pointing

class Telescope:

    _id = 0
    def __init__(self, x, y, z, focal, camera_radius, toDeg=False):
        
        Telescope._id += 1
        self.id = Telescope._id
        self.x = x.to(u.m)
        self.y = y.to(u.m)
        self.z = z.to(u.m)
        self.focal = focal.to(u.m)
        self.camera_radius = camera_radius.to(u.m)
        self.alt = u.Quantity(0, u.rad)
        self.az = u.Quantity(0, u.rad)
        
        self.table = self.make_table(toDeg=toDeg)
        
    def point_to_altaz(self, alt, az):
        self.alt = alt.to(u.rad)
        self.az = az.to(u.rad)
        if self.az < 0:
            self.az += 2*np.pi*u.rad
        self.update_table(["alt", "az", "zn", "pointing"])

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

    def point_to_object(self, object):
        """
        Point to object.

        Parameters
        ----------
        object: numpy.array([x, y, z])
        """
        GT = np.sqrt(((self.position - object) ** 2).sum())
        alt_tel = np.arcsin((-self.z.value + object[2]) / GT)
        az_tel = np.arctan2((- self.y.value + object[1]), (- self.x.value + object[0]))
        self.point_to_altaz(alt_tel * u.rad, az_tel * u.rad)
        self.update_table(["alt", "az", "zn", "pointing"])

    @property
    def pointing_vector(self):
        # return pointing.alt_az_to_vector(self.alt, self.az)
        return np.array([np.cos(self.alt.to(u.rad))*np.cos(self.az.to(u.rad)),
                         -np.cos(self.alt.to(u.rad))*np.sin(self.az.to(u.rad)),
                         np.sin(self.alt.to(u.rad))])

    def make_table(self, toDeg=False):
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
        
        if toDeg:
            table = utils.deg2rad(table, toDeg=toDeg)
            return table
        else:
            return table

    def update_table(self, vals):
        for val in vals:
            if val == "pointing":
                new_val = self.pointing_vector
                for i, v in enumerate(["p_x", "p_y", "p_z"]):
                    self.table[v] = new_val[i]
                    self.table[v].info.format = '{:.3f}'
                
            else:
                self.table[val] = getattr(self, val)
                self.table[val].info.format = '{:.3f}'

    def convert_unit(self, toDeg=True):
        self.table = utils.deg2rad(self.table, toDeg)
        

class Array:

    def __init__(self, telescope_list, frame=None, **kwargs):
        
        self.telescopes = telescope_list
        self._div = 0
        self._pointing = {"az":0*u.deg, "alt":0*u.deg}

        if frame == None:
            self._frame = utils.CTA_INFO(verbose=False, **kwargs)

        self.__create_table__()

    def __create_table__(self):
        table = []
        for tel in self.telescopes:
            table.append(tel.table)
        
        self.table = tab.vstack(table)

        self.table.add_column(self.dist2tel, name="d_tel")
        self.table["d_tel"].unit = u.m
        self.table["d_tel"].info.format = "{:.2f}"

    def convert_unit(self, toDeg=True, showTable=True):
        self.table = utils.deg2rad(self.table, toDeg)
        if showTable:
            return self.table

    @property
    def frame(self):
        return self._frame
    
    @property
    def positions_array(self):
        """
        all telescopes positions as an array

        Returns
        -------
        np.array
        """
        return np.array([tel.position for tel in self.telescopes])

    @property
    def pointing_vectors(self):
        """
        all telescopes pointing vectors as an array

        Returns
        -------
        np.array
        """
        return np.array([tel.pointing_vector for tel in self.telescopes])
        

    @property
    def positions_table(self):
        """
        all telescopes positions as a table

        Returns
        -------
        astropy.table
        """
        return self.table["x", "y", "z"]
        
    
    @property        
    def pointing_table(self):
        """
        all telescopes pointing vectors as a table

        Returns
        -------
        astropy.table
        """
        return self.table["p_x", "p_y", "p_z"]
    
    
    def pointing_coord(self, icrs=True):
        return pointing.pointing_coord(self.table, self.frame, icrs=icrs)
    
    @property
    def barycenter(self):
        return np.array(utils.calc_mean(self.table, ["x", "y", "z"]))

    @property
    def dist2tel(self):
        dist = np.zeros(self.table.__len__())
        for i, axis in enumerate(["x", "y", "z"]):
            dist += (self.table[axis] - self.barycenter[i])**2.
        dist = np.sqrt(dist)
        return dist
        #return np.sqrt(np.sum((self.positions_array - self.barycenter)**2, axis=1))

    @property
    def div(self):
        return self._div

    @property
    def pointing(self):
        return self._pointing
    
    def divergent_pointing(self, div, alt_mean, az_mean):
        """
        Divergent pointing given a parameter div.
        Update pointing of all telescopes of the array

        Parameters
        ----------
        div: float between 0 and 1
        alt_mean: `astropy.Quantity`
            mean alt pointing
        az_mean: `astropy.Quantity`
            mean az pointing
        """
        if type(alt_mean) != u.Quantity:
            print("[Warning] The unit of the input alt is assumed to be 'deg'")
            alt_mean = alt_mean*u.deg

        if type(az_mean) != u.Quantity:
            print("[Warning] The unit of the input az is assumed to be 'deg'")
            az_mean = az_mean*u.deg
            if az_mean.value < 0:
                az_mean += 360*u.deg

        self._div = div
        self._pointing["az"] = az_mean
        self._pointing["alt"] = alt_mean

        if div > 1 or div < 0:
            print("[Error] The div value should be between 0 and 1.")
        elif div==0:
            for tel in self.telescopes:
                tel.point_to_altaz(alt_mean, az_mean)
        else:
            G = pointing.pointG_position(self.barycenter, div, alt_mean, az_mean)
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
        

        if skymap:
            if group:
                tel_group, color = utils.group_table(self.table, group)
                
                for i, [table, c] in enumerate(zip(tel_group.groups, color)):
                
                    radec = pointing.pointing_coord(table, self.frame, icrs=True)
                    ax = visual.display_skymap(radec, self.frame,  
                                        label="group_{}".format(i), ax=ax, color=c,)
            else:
                radec = pointing.pointing_coord(self.table, self.frame, icrs=True)
                ax = visual.display_skymap(radec, self.frame, ax=ax)
            return ax
        else:
            proj = []
            for axis in ["x", "y", "z"]:
                if axis in projection:
                    proj.append(axis)

        if len(proj) == 1:
            ax = visual.display_1d(self.table, proj, group=group, ax=ax, **kwargs)
        elif len(proj) == 2:
            ax = visual.display_2d(self.table, proj, group=group, ax=ax, **kwargs)
        else:
            ax = visual.display_3d(self.table, proj, group=group, ax=ax, **kwargs)
        
        return ax

    def skymap_polar(self, group, fig=None):
        return visual.skymap_polar(self, group=group, fig=fig)

    def multiplicity_plot(self, fig=None):
        return visual.multiplicity_plot(self, fig=fig)

    def export_cfg(self, filename=None):
        
        #create  simtel cfg file
        if filename==None:
            filename = 'CTA-ULTRA6-LaPalma-div{}-az{}-alt{}.cfg'.format(
                str(self.div).replace(".", "_"), 
                str(self.pointing["az"].value).replace(".", "_"), 
                str(self.pointing["alt"].value).replace(".", "_"))
        
        with open(filename, 'w') as f:
            f.write('#ifndef TELESCOPE')
            f.write('# define TELESCOPE 0')
            f.write('#endif\n')
            f.write('#if TELESCOPE == 0\n""")')
            f.write('   TELESCOPE_THETA={:.2f} \n'.format(90 - self.pointing["alt"].value))
            f.write('   TELESCOPE_PHI={:.2f} \n'.format(self.pointing["az"].value))
            f.write("""\n% Global and default configuration for things missing in telescope-specific config.
                # include <CTA-ULTRA6-LST.cfg>\n""")
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
            f.write('Error Invalid telescope for CTA-ULTRA6 La Palma configuration.\n')
            f.write('#endif\n')
            f.write('trigger_telescopes = 2 % In contrast to Prod-3 South we apply loose stereo trigger immediately')
            f.write('array_trigger = array_trigger_ultra6_diver-test.dat')
        

