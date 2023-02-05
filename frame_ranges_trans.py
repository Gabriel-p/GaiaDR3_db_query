
import pandas as pd
import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord


def radec2lonlat(ra, dec):
    gc = SkyCoord(ra=ra * u.degree, dec=dec * u.degree)
    lb = gc.transform_to('galactic')
    lon, lat = lb.l.value, lb.b.value
    return lon, lat


df = pd.read_csv("frame_ranges.txt")
ra_min, ra_max, dec_min, dec_max = df['ra_min'], df['ra_max'], df['dec_min'],\
    df['dec_max']

# Galactic quadrant
lon1, lat1 = radec2lonlat(ra_min, dec_min)
lon2, lat2 = radec2lonlat(ra_min, dec_max)
lon3, lat3 = radec2lonlat(ra_max, dec_min)
lon4, lat4 = radec2lonlat(ra_max, dec_max)

lmin = np.min([lon1, lon2, lon3, lon4], axis=0)
lmax = np.max([lon1, lon2, lon3, lon4], axis=0)
bmin = np.min([lat1, lat2, lat3, lat4], axis=0)
bmax = np.max([lat1, lat2, lat3, lat4], axis=0)

df['ra_min'] = np.round(lmin, 3)
df['ra_max'] = np.round(lmax, 3)
df['dec_min'] = np.round(bmin, 3)
df['dec_max'] = np.round(bmax, 3)

df.rename(columns={
    'ra_min': 'l_min', 'ra_max': 'l_max', 'dec_min': 'b_min',
    'dec_max': 'b_max'}, inplace=True)

df.to_csv("frame_ranges_galac.txt", index=False)
