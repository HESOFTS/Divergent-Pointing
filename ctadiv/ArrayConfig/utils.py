import astropy.units as u

def deg2rad(table, deg=False):
    if deg:
        table["alt"]=table["alt"].to(u.deg)
        table["az"]=table["az"].to(u.deg)
        table["zn"]=table["zn"].to(u.deg)
        table["fov"]=table["fov"].to(u.deg**2)
    else:
        table["alt"]=table["alt"].to(u.rad)
        table["az"]=table["az"].to(u.rad)
        table["zn"]=table["zn"].to(u.rad)
        table["fov"]=table["fov"].to(u.rad**2)
    return table