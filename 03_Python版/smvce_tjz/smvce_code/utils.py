"""
Utility functions for SM-VCE toolkit.
Converted from MATLAB: splitstr.m, lonlat2sub.m, sub2lonlat.m, coorlim.m, datalim.m, coor2dempar.m
"""
import numpy as np


def splitstr(line):
    """Split a string by commas, colons, spaces, and tabs."""
    separators = set([',', ':', ' ', '\t'])
    parts = []
    current = []
    for ch in line:
        if ch in separators:
            if current:
                parts.append(''.join(current))
                current = []
        else:
            current.append(ch)
    if current:
        parts.append(''.join(current))
    return parts


def lonlat2sub(lonlat, coor):
    """
    Get the sub index of lonlat based on the coordinate of coor.

    Parameters
    ----------
    lonlat : ndarray, shape (N, 2) - [lon, lat]
    coor : dict with keys corner_lon, corner_lat, post_lon, post_lat

    Returns
    -------
    subs : ndarray, shape (N, 2) - [row, col] (0-based indexing)
    """
    lonlat = np.atleast_2d(lonlat)
    rows = np.round((lonlat[:, 1] - coor['corner_lat']) / coor['post_lat']).astype(int)
    cols = np.round((lonlat[:, 0] - coor['corner_lon']) / coor['post_lon']).astype(int)
    return np.column_stack([rows, cols])


def sub2lonlat(subs, coor):
    """
    Convert sub indices to lon/lat coordinates.

    Parameters
    ----------
    subs : ndarray, shape (N, 2) - [row, col] (0-based indexing)
    coor : dict

    Returns
    -------
    lonlat : ndarray, shape (N, 2) - [lon, lat]
    """
    subs = np.atleast_2d(subs)
    lons = subs[:, 1] * coor['post_lon'] + coor['corner_lon']
    lats = subs[:, 0] * coor['post_lat'] + coor['corner_lat']
    return np.column_stack([lons, lats])


def coorlim(coor):
    """Calculate coordinate limits."""
    return [
        coor['corner_lon'],
        coor['corner_lon'] + coor['post_lon'] * (coor['width'] - 1),
        coor['corner_lat'] + coor['post_lat'] * (coor['nlines'] - 1),
        coor['corner_lat']
    ]


def datalim(data, ord_val=3):
    """Calculate appropriate data limits for plotting."""
    d = data.copy()
    d[np.isnan(d)] = 0
    d[np.abs(d) == np.inf] = 0
    a = d[d != 0]
    if len(a) == 0:
        return [-0.5, 0.5]
    m = np.nanmean(a)
    sd = np.nanstd(a)
    return [m - ord_val * sd, m + ord_val * sd]


def coor2dempar(coor, pout):
    """Write coordinate information to a DEM parameter file."""
    with open(pout, 'w') as f:
        f.write('Gamma DIFF&GEO DEM/MAP parameter file\n')
        f.write('title: DEM\n')
        f.write('DEM_projection:     EQA\n')
        f.write('data_format:        REAL*4\n')
        f.write('DEM_hgt_offset:          0.00000\n')
        f.write('DEM_scale:               1.00000\n')
        f.write(f'width:                {coor["width"]}\n')
        f.write(f'nlines:               {coor["nlines"]}\n')
        f.write(f'corner_lat:     {coor["corner_lat"]:.7f}  decimal degrees\n')
        f.write(f'corner_lon:   {coor["corner_lon"]:.7f}  decimal degrees\n')
        f.write(f'post_lat:   {coor["post_lat"]:.7f}  decimal degrees\n')
        f.write(f'post_lon:    {coor["post_lon"]:.7f}  decimal degrees\n')
        f.write('\n')
        f.write('ellipsoid_name: WGS 84\n')
        f.write('ellipsoid_ra:        6378137.000   m\n')
        f.write('ellipsoid_reciprocal_flattening:  298.2572236\n')
        f.write('\n')
        f.write('datum_name: WGS 1984\n')
        f.write('datum_shift_dx:              0.000   m\n')
        f.write('datum_shift_dy:              0.000   m\n')
        f.write('datum_shift_dz:              0.000   m\n')
        f.write('datum_scale_m:         0.00000e+00\n')
        f.write('datum_rotation_alpha:  0.00000e+00   arc-sec\n')
        f.write('datum_rotation_beta:   0.00000e+00   arc-sec\n')
        f.write('datum_rotation_gamma:  0.00000e+00   arc-sec\n')
        f.write('datum_country_list: Global Definition, WGS84, World\n')
