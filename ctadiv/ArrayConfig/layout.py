from ..const import CONFIG_DIR
from .telescope import Telescope, Array
import astropy.units as u
import numpy as np

from . import utils

def LoadConfig(file, tel_id=-1, radius="degrees", frame=None, **kwargs):
    
    """
    Load the telescope configuration file
    
    Parameters
    ----------
    file: str 
          the name of file
    tel_id: int 
            If you want to load only a single telescope,
            you can set this parameter (defalut: load all)
    radius: str
            Define the unit of camera radius
            either 'meters' or 'degrees'
        
    Returns
    -------
    class.Array
    """

    with open(file, "r") as f:
        
        tels = []
        for line in f.readlines():
            line = np.asarray(line.split()).astype("float")

            if (radius!="meters"):
                line[4] = utils.convert_radius(line[4]*u.deg, line[3], toDeg=False)
            
            coord = [x*u.m for x in line]

            tel = Telescope(coord[0],coord[1],coord[2],coord[3],coord[4])
            tels.append(tel)

    if tel_id == -1:
        return Array(tels, frame=frame, **kwargs)
    else:
        return tels[tel_id-1]