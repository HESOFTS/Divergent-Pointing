from astropy.table import Table
import astropy.units as u

import numpy as np
from ..const import SCRIPT_DIR

# Load FoV
class loadFov():
    _listOfFoV = ["05", "2", "4"]
    
    def __init__(self):
        self.__loadFov__()

        for mode in self.listOfFoV:
            setattr(self, "FoV_"+str(mode), self.__selectFoV__(mode))

    def __repr__(self):
        return self.table.__repr__()

    @property
    def table(self):
        return self._table
    
    @property
    def listOfFoV(self):
        return self._listOfFoV
    
    def __loadFov__(self):
        tels_dict = []
        for mode in self.listOfFoV:
            with open(SCRIPT_DIR+"/config/CTA-ULTRA6-LaPalma-divergent_{}_180.cfg".format(mode)) as div:
                text = div.read()
                text = text.split("#")[1:]

                for line in text:
                    line_list = line.split("\n")
                    tels_dict.append([mode, float(line_list[1].split("=")[1]),  float(line_list[2].split("=")[1])])

        self._table = Table(np.asarray(tels_dict), names=["FoV", "Theta", "Phi"], units=["", u.deg, u.deg])

    def __selectFoV__(self, mode):
        return self.table[self.table["FoV"]==mode]