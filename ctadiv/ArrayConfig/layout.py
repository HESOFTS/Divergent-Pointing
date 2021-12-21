from ..const import CONFIG_DIR
from .telescope import Telescope, Array
import astropy.units as u

def LoadLayout(file="layout-3AL4M15-5.txt", tel_id=-1):

    with open(CONFIG_DIR+file, "r") as f:
        
        tels = []
        for line in f:  
            #split the string on whitespace, return a list of numbers as strings
            coord_str = line.split()
            coord_str[0], coord_str[1], coord_str[2] = float(coord_str[0]), float(coord_str[1]), float(coord_str[2]) 
            coord = [x*u.m for x in coord_str]
            tel = Telescope(coord[0],coord[1],coord[2],coord[3],coord[4])
            tels.append(tel)

    if tel_id == -1:
        return Array(tels)
    else:
        return tels[tel_id-1]