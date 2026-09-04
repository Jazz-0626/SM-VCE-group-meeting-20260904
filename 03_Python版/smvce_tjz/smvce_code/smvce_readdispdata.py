"""
Read displacement data for SM-VCE processing.
Converted from MATLAB: SMVCE_readdispdata.m
"""
import os
import numpy as np
import rasterio


def smvce_readdispdata(pdata=None):
    """
    Read the prepared SAR displacement data.

    Parameters
    ----------
    pdata : str, optional - path to data folder (default: ./SMVCE_DATA)

    Returns
    -------
    data : ndarray, shape (row, col, N) - displacement measurements
    inc : ndarray, shape (row, col, N) - incidence angles
    azi : ndarray, shape (row, col, N) - azimuth angles
    losazienu : ndarray, shape (N,) - observation geometry flag
    leftorright : ndarray, shape (N,) - looking mode
    coor : dict - coordinate information
    dem : ndarray, shape (row, col) - DEM data
    mask : ndarray, shape (row, col) - mask
    fault : list or 0 - fault traces
    datainfo : list - data information
    """
    if pdata is None:
        pdata = os.path.join(os.getcwd(), 'SMVCE_DATA')

    pdata_information = os.path.join(pdata, 'data_information')

    # Read data_information file
    datainfo = []
    with open(pdata_information, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = _splitstr(line)
            datainfo.append(parts)

    N = len(datainfo)
    losazienu = np.zeros(N, dtype=int)
    leftorright = np.zeros(N, dtype=int)

    data_list = []
    inc_list = []
    azi_list = []
    coor = None

    for i in range(N):
        losazienu[i] = int(datainfo[i][1])
        leftorright[i] = int(datainfo[i][4])

        ptifi = os.path.join(pdata, datainfo[i][0])
        datai, transform, crs = _read_geotiff(ptifi)
        data_list.append(datai)

        if i == 0:
            row, col = datai.shape
            coor = {
                'corner_lon': transform.c,
                'corner_lat': transform.f,
                'post_lon': transform.a,
                'post_lat': transform.e,
                'nlines': row,
                'width': col,
                'crs': crs.to_string() if crs is not None else None
            }

        row_i, col_i = datai.shape

        # Read or set incidence angle
        try:
            inc_val = float(datainfo[i][2])
            inci = np.full((row_i, col_i), inc_val)
        except (ValueError, IndexError):
            inc_file = os.path.join(pdata, datainfo[i][2])
            inci, _, _ = _read_geotiff(inc_file)

        inc_list.append(inci)

        # Read or set azimuth angle
        try:
            azi_val = float(datainfo[i][3])
            azii = np.full((row_i, col_i), azi_val)
        except (ValueError, IndexError):
            azi_file = os.path.join(pdata, datainfo[i][3])
            azii, _, _ = _read_geotiff(azi_file)

        azi_list.append(azii)

    data = np.stack(data_list, axis=2)
    inc = np.stack(inc_list, axis=2)
    azi = np.stack(azi_list, axis=2)
    row, col, _ = data.shape

    # Read DEM
    dem_path = os.path.join(pdata, 'dem.tif')
    if os.path.exists(dem_path):
        dem, _, _ = _read_geotiff(dem_path)
    else:
        dem = np.random.randn(row, col)

    # Read mask
    mask_path = os.path.join(pdata, 'mask.tif')
    if os.path.exists(mask_path):
        mask, _, _ = _read_geotiff(mask_path)
    else:
        mask = np.ones((row, col))

    # Read fault
    fault_path = os.path.join(pdata, 'fault.xy')
    if os.path.exists(fault_path):
        fault = _readfault(fault_path)
    else:
        fault = 0

    return data, inc, azi, losazienu, leftorright, coor, dem, mask, fault, datainfo


def _read_geotiff(filepath):
    """Read a GeoTIFF file and return data, transform, and CRS."""
    with rasterio.open(filepath) as src:
        data = src.read(1).astype(np.float64)
        transform = src.transform
        crs = src.crs
    return data, transform, crs


def _readfault(pfault):
    """Read fault trace data from a file."""
    fault = []
    current_coor = None

    with open(pfault, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line == ' ':
                continue
            if '>' in line:
                if current_coor is not None and len(current_coor) > 0:
                    fault.append(np.array(current_coor))
                current_coor = []
                continue
            parts = line.split(',')
            if len(parts) >= 2:
                try:
                    lon = float(parts[0])
                    lat = float(parts[1])
                    if current_coor is None:
                        current_coor = []
                    current_coor.append([lon, lat])
                except ValueError:
                    continue

    if current_coor is not None and len(current_coor) > 0:
        fault.append(np.array(current_coor))

    return fault if fault else 0


def _splitstr(line):
    """Split string by separators (comma, colon, space, tab)."""
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
