
import pathlib
import gzip
import warnings
import pandas as pd
import astropy.units as u
# from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
import numpy as np
# from uncertainties import ufloat
# from uncertainties import unumpy as unp
import matplotlib.pyplot as plt


# ******** IMPORTANT ********
#
# Clusters that wrap around the edges of the (ra, dec) coordinates are not
# still properly process; e.g.: Blanco 1
#
# ******** IMPORTANT ********


# # Gaia EDR3 zero points. Sigmas are already squared here.
Zp_G, sigma_ZG_2 = 25.6873668671, 0.00000759
Zp_BP, sigma_ZBP_2 = 25.3385422158, 0.000007785
Zp_RP, sigma_ZRP_2 = 24.7478955012, 0.00001428

# Path to the database
db_path = '/media/gabriel/backup/gabriel/GaiaDR3/datafiles/'

# Input file
input_f = "clusters.ini"

# File that contains the regions delimited by each frame
txt_f = 'frame_ranges.txt'

#
# Input parameters
frame = 'galactic'  # 'equatorial'
# Mag (flux) filter
max_mag = 19
# data['Gmag'] = Zp_G + -2.5 * np.log10(I_G)
min_G_flux = 10**((max_mag - Zp_G) / (-2.5))
# Maximum size of box to query (in degrees)
max_box = 50


def main():
    """
    box_s is assumed to be in Equatorial coordinates
    """
    print(f"Using {frame} frame")

    clusters = pd.read_csv(input_f)

    # Read data about the regions occupied by each frame
    fdata = pd.read_csv(txt_f)

    for i, cluster in clusters.iterrows():
        parseCluster(fdata, cluster)


def parseCluster(fdata, cluster):
    """
    """
    cl_name = cluster['Name']
    print(f"\nQuerying cluster: {cl_name}...")

    # Cluster coordinates: center and box size
    c_ra, c_dec = cluster['ra'], cluster['dec']

    # lon_lat_c = radec2lonlat(c_ra, c_dec)
    # if lon_lat_c[0] + cluster['box'] > 360 or lon_lat_c[0] - cluster['box'] < 0:
    #     print(cl_name, lon_lat_c, cluster['box'])
    # return

    # Size of box to query (in degrees)
    box_s_eq = cluster['box']

    data_in_files, xmin_cl, xmax_cl, ymin_cl, ymax_cl = findFrames(
        c_ra, c_dec, box_s_eq, fdata)

    if len(data_in_files) == 0:
        return

    plx_d = float(cluster['plx_d'])
    if plx_d > 10:
        plx_p = .5
    elif plx_d <= 10 and plx_d > 5:
        plx_p = .3
    else:
        plx_p = .2
    plx_min = plx_d - plx_d * plx_p

    all_frames = query(
        cluster, data_in_files, xmin_cl, xmax_cl, ymin_cl, ymax_cl,
        plx_min)

    print("Adding uncertainties")
    # all_frames = xyCoords(all_frames)
    all_frames = uncertMags(all_frames)
    all_frames = all_frames.drop(columns=[
        'FG', 'e_FG', 'FBP', 'e_FBP', 'FRP', 'e_FRP'])

    # msk = all_frames['Gmag'] < max_mag
    # all_frames = all_frames[msk]
    # print(f"{len(all_frames)} stars after magnitude cut <{max_mag}")

    print("Save to file")
    all_frames.to_csv('./out/' + cl_name + '.csv', index=False)


def findFrames(c_ra, c_dec, box_s_eq, fdata):
    """
    """
    # These are the points that determine the range of *all* the frames
    ra_min, ra_max = fdata['ra_min'], fdata['ra_max']
    dec_min, dec_max = fdata['dec_min'], fdata['dec_max']

    box_s_eq = min(max_box, box_s_eq)
    print("({:.3f}, {:.3f}); Box size: {:.2f} deg".format(
        c_ra, c_dec, box_s_eq))

    if frame == 'galactic':
        box_s_eq = np.sqrt(2) * box_s_eq

    # Correct size in RA
    box_s_x = box_s_eq / np.cos(np.deg2rad(c_dec))

    xl, yl = box_s_x * .5, box_s_eq * .5

    # Limits of the cluster's region in Equatorial
    xmin_cl, xmax_cl = c_ra - xl, c_ra + xl
    ymin_cl, ymax_cl = c_dec - yl, c_dec + yl

    # Identify which frames contain the cluster region
    l2 = (xmin_cl, ymax_cl)
    r2 = (xmax_cl, ymin_cl)
    frame_intersec = []
    for i, xmin_fr_i in enumerate(ra_min):
        # Top left
        l1 = (xmin_fr_i, dec_max[i])
        # Bottom right
        r1 = (ra_max[i], dec_min[i])
        frame_intersec.append(doOverlap(l1, r1, l2, r2))

    frame_intersec = np.array(frame_intersec)

    data_in_files = list(fdata[frame_intersec]['filename'])
    print(f"Cluster is present in {len(data_in_files)} frames")

    # plt.scatter(xmin_cl, ymin_cl, c='k')
    # plt.scatter(xmin_cl, ymax_cl, c='k')
    # plt.scatter(xmax_cl, ymin_cl, c='k')
    # plt.scatter(xmax_cl, ymax_cl, c='k')
    # from matplotlib.pyplot import cm
    # color = iter(cm.rainbow(np.linspace(0, 1, len(data_in_files))))
    # for i, flag in enumerate(frame_intersec):
    #     if flag:
    #         c = next(color)
    #         plt.scatter(xmin_fr[i], ymin_fr[i], color=c, marker='x')
    #         plt.scatter(xmin_fr[i], ymax_fr[i], color=c, marker='x')
    #         plt.scatter(xmax_fr[i], ymax_fr[i], color=c, marker='x')
    #         plt.scatter(xmax_fr[i], ymin_fr[i], color=c, marker='x')
    # plt.show()

    return data_in_files, xmin_cl, xmax_cl, ymin_cl, ymax_cl


def doOverlap(l1, r1, l2, r2):
    """
    Given two rectangles, find if the given two rectangles overlap or not.
    Note that a rectangle can be represented by two coordinates, top left and
    bottom right.
    l1: Top Left coordinate of first rectangle.
    r1: Bottom Right coordinate of first rectangle.
    l2: Top Left coordinate of second rectangle.
    r2: Bottom Right coordinate of second rectangle.

    Source: https://www.geeksforgeeks.org/find-two-rectangles-overlap/
    """
    min_x1, max_y1 = l1
    max_x1, min_y1 = r1
    min_x2, max_y2 = l2
    max_x2, min_y2 = r2
    # If one rectangle is on left side of other
    if (min_x1 > max_x2 or min_x2 > max_x1):
        return False
    # If one rectangle is above other
    if(min_y1 > max_y2 or min_y2 > max_y1):
        return False
    return True


def query(
    cluster, data_in_files, xmin_cl, xmax_cl, ymin_cl, ymax_cl,
        plx_min):
    """
    """
    all_frames = []
    for i, file in enumerate(data_in_files):
        with gzip.open(db_path + file) as f:
            data = pd.read_csv(f, index_col=False)

            mx = (data['ra'] >= xmin_cl) & (data['ra'] <= xmax_cl)
            my = (data['dec'] >= ymin_cl) & (data['dec'] <= ymax_cl)
            m_plx = data['parallax'] > plx_min
            m_gmag = data['phot_g_mean_flux'] > min_G_flux
            msk = (mx & my & m_plx & m_gmag)

            print(f"{i+1}, {file} contains {msk.sum()} cluster stars")
            if msk.sum() == 0:
                continue

            # data['lon'], data['lat'] = radec2lonlat(
            #     data['ra'].values, data['dec'].values)

            all_frames.append(data[msk])
    all_frames = pd.concat(all_frames)

    if frame == 'galactic':
        c_ra, c_dec = cluster['ra'], cluster['dec']
        box_s_h = cluster['box'] * .5
        gal_cent = radec2lonlat(c_ra, c_dec)

        if all_frames['l'].max() - all_frames['l'].min() > 180:
            warnings.warn("Frame wraps around 360 in longitude. Fixing..")
            breakpoint()

            lon = all_frames['l'].values
            if gal_cent[0] > 180:
                msk = lon < 180
                lon[msk] += 360
            else:
                msk = lon > 180
                lon[msk] -= 360
            all_frames['l'] = lon

        xmin_cl, xmax_cl = gal_cent[0] - box_s_h, gal_cent[0] + box_s_h
        ymin_cl, ymax_cl = gal_cent[1] - box_s_h, gal_cent[1] + box_s_h
        mx = (all_frames['l'] >= xmin_cl) & (all_frames['l'] <= xmax_cl)
        my = (all_frames['b'] >= ymin_cl) & (all_frames['b'] <= ymax_cl)
        msk = (mx & my)
        all_frames = all_frames[msk]

    all_frames = all_frames.rename(columns={
        'source_id': 'Source', 'ra': 'RA_ICRS', 'dec': 'DE_ICRS',
        'parallax': 'Plx', 'parallax_error': 'e_Plx',
        'pmra': 'pmRA', 'pmra_error': 'e_pmRA', 'b': 'GLAT',
        'pmdec': 'pmDE', 'pmdec_error': 'e_pmDE', 'l': 'GLON',
        'phot_g_mean_flux': 'FG', 'phot_g_mean_flux_error': 'e_FG',
        'phot_bp_mean_flux': 'FBP', 'phot_bp_mean_flux_error': 'e_FBP',
        'phot_rp_mean_flux': 'FRP', 'phot_rp_mean_flux_error': 'e_FRP',
        'radial_velocity': 'RV', 'radial_velocity_error': 'e_RV'}
    )

    print(f"{len(all_frames)} stars queried")
    return all_frames


def xyCoords(data):
    """
    Replicate Vizier's '_x, _y' columns
    """
    ra, de = data['RA_ICRS'], data['DE_ICRS']
    # Center of the frame
    ra_0, de_0 = .5 * (np.max(ra) + np.min(ra)), .5 * (np.max(de) + np.min(de))
    coord1 = SkyCoord(ra_0, de_0, frame='icrs', unit='deg')
    coord2 = SkyCoord(ra.values, de.values, frame='icrs', unit='deg')
    r = coord1.separation(coord2)
    theta = coord1.position_angle(coord2)
    _x = r * np.cos(theta)
    _y = r * np.sin(theta)

    data['_x'], data['_y'] = _x.value, _y.value
    # data = data.assign(_x=_x)
    # data = data.assign(_y=_y)

    return data


def uncertMags(data):
    """
    # Gaia EDR3 zero points:

    https://www.cosmos.esa.int/web/gaia/edr3-passbands

    """
    I_G, e_IG = data['FG'].values, data['e_FG'].values
    I_BP, e_IBP = data['FBP'].values, data['e_FBP'].values
    I_RP, e_IRP = data['FRP'].values, data['e_FRP'].values

    data['Gmag'] = Zp_G + -2.5 * np.log10(I_G)
    BPmag = Zp_BP + -2.5 * np.log10(I_BP)
    RPmag = Zp_RP + -2.5 * np.log10(I_RP)
    data['BP-RP'] = BPmag - RPmag

    e_G = np.sqrt(sigma_ZG_2 + 1.179 * (e_IG / I_G)**2)
    data['e_Gmag'] = e_G
    e_BP = np.sqrt(sigma_ZBP_2 + 1.179 * (e_IBP / I_BP)**2)
    e_RP = np.sqrt(sigma_ZRP_2 + 1.179 * (e_IRP / I_RP)**2)
    data['e_BP-RP'] = np.sqrt(e_BP**2 + e_RP**2)

    return data


def radec2lonlat(ra, dec):
    gc = SkyCoord(ra=ra * u.degree, dec=dec * u.degree)
    lb = gc.transform_to('galactic')
    lon, lat = lb.l.value, lb.b.value
    return lon, lat


if __name__ == '__main__':
    # Create output dir if it does not exist
    path = pathlib.Path('./out')
    path.mkdir(parents=True, exist_ok=True)
    main()
