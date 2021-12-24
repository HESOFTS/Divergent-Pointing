from ..const import CONFIG_DIR
from .telescope import Telescope, Array
import astropy.units as u
import numpy as np

from . import utils

def LoadConfig(file, tel_id=-1, radius="deg"):

    with open(file, "r") as f:
        
        tels = []
        for line in f.readlines():
            line = np.asarray(line.split()).astype("float")

            if radius=="deg":
                line[4] = utils.convert_radius(line[4]*u.deg, line[3], toDeg=False)
            
            coord = [x*u.m for x in line]

            tel = Telescope(coord[0],coord[1],coord[2],coord[3],coord[4])
            tels.append(tel)

    if tel_id == -1:
        return Array(tels)
    else:
        return tels[tel_id-1]