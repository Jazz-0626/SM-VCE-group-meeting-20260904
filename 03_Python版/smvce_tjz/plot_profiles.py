"""Combined profile figure: E-W map (left) + N equally-spaced ENU profiles (right).

make_profile_figure(result_dir) reads enu_{east,north,up}.tif and draws:
  - left : E-W displacement map with the profile lines overlaid
  - right: N stacked subplots, each the SM-VCE E-W / N-S / Vertical cross-section
           (W->E) along that profile.

Profiles span LAT0..LAT1 (near-field), equally spaced.
"""
import os, sys
os.environ.pop('PROJ_DATA', None); os.environ.pop('PROJ_LIB', None); os.environ.pop('GDAL_DATA', None)
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import TwoSlopeNorm
import rasterio
from scipy.ndimage import binary_erosion

LAT0, LAT1, NPROF = 28.65, 28.90, 10     # near-field profile band (former P3..P8)
FAULT_LON = 87.55                         # approx. fault trace longitude (marker)

def _rd(path):
    with rasterio.open(path) as s:
        a = s.read(1).astype(float)
        if s.nodata is not None:
            a[a == s.nodata] = np.nan
        tr = s.transform
    return a, tr

def _sample(field, tr, lons, lat):
    r, c = field.shape
    x0, y0 = tr.c, tr.f
    dx, dy = abs(tr.a), abs(tr.e)
    out = np.full(len(lons), np.nan)
    ri = int(round((y0 - lat) / dy))
    if not (0 <= ri < r):
        return out
    for k, lon in enumerate(lons):
        ci = int(round((lon - x0) / dx))
        if 0 <= ci < c:
            v = field[ri, ci]
            if np.isfinite(v):
                out[k] = v
    return out

def make_profile_figure(result_dir):
    g = os.path.join(result_dir, 'geotiff')
    e, tr = _rd(os.path.join(g, 'enu_east.tif'))
    n, _ = _rd(os.path.join(g, 'enu_north.tif'))
    u, _ = _rd(os.path.join(g, 'enu_up.tif'))
    H, W = e.shape
    lonW, lonE = tr.c, tr.c + W * tr.a
    latN, latS = tr.f, tr.f + H * tr.e

    # cleaned E-W for the map background
    solved = np.isfinite(e) & np.isfinite(n) & np.isfinite(u)
    core = binary_erosion(solved, iterations=10)
    ew_map = np.where(core, e, np.nan)
    ew_map[np.abs(ew_map) > 1.5] = np.nan

    lats = np.linspace(LAT0, LAT1, NPROF)
    colors = plt.cm.turbo(np.linspace(0.05, 0.95, NPROF))
    lons = np.linspace(lonW, lonE, 400)
    dist = (lons - lonW) * 111.32 * np.cos(np.deg2rad((LAT0 + LAT1) / 2))
    fault_d = (FAULT_LON - lonW) * 111.32 * np.cos(np.deg2rad((LAT0 + LAT1) / 2))

    fig = plt.figure(figsize=(15, 13))
    gs = gridspec.GridSpec(NPROF, 2, width_ratios=[1.0, 1.5], wspace=0.22, hspace=0.18)

    # ---- left: E-W map + profile lines ----
    axm = fig.add_subplot(gs[:, 0])
    norm = TwoSlopeNorm(vmin=-1.0, vcenter=0, vmax=1.0)
    axm.imshow(ew_map, extent=[lonW, lonE, latS, latN], origin='upper',
               cmap='RdBu_r', norm=norm, aspect='auto', interpolation='bilinear')
    axm.axvline(FAULT_LON, color='k', ls='--', lw=0.8, alpha=0.5)
    for k, lat in enumerate(lats):
        axm.axhline(lat, color=colors[k], lw=1.8)
        axm.text(lonW + 0.01, lat + 0.006, f'P{k+1}', color=colors[k],
                 fontsize=9, fontweight='bold', va='bottom')
    axm.plot(87.45, 28.50, 'k*', ms=9)
    axm.set_xlim(lonW, lonE); axm.set_ylim(latS, latN)
    axm.set_xlabel('Longitude (°E)'); axm.set_ylabel('Latitude (°N)')
    axm.set_title('E-W Displacement + 10 Profiles', fontweight='bold')

    # ---- right: stacked ENU profiles ----
    for k, lat in enumerate(lats):
        axp = fig.add_subplot(gs[k, 1])
        ep = _sample(e, tr, lons, lat); npf = _sample(n, tr, lons, lat); up = _sample(u, tr, lons, lat)
        axp.axhline(0, color='gray', ls=':', lw=0.5)
        axp.axvline(fault_d, color='red', ls=':', lw=1.0, alpha=0.6)
        axp.plot(dist, ep, '-', color='black', lw=1.6, label='E-W')
        axp.plot(dist, up, '-', color='#c0392b', lw=1.4, label='Up')
        axp.plot(dist, npf, '-', color='#6a3d9a', lw=1.2, alpha=0.85, label='N-S')
        axp.set_ylim(-2.6, 1.2)
        axp.set_ylabel(f'P{k+1}', color=colors[k], fontweight='bold', fontsize=9, rotation=0, labelpad=14, va='center')
        axp.tick_params(labelsize=7)
        axp.grid(True, alpha=0.15)
        axp.text(0.99, 0.92, f'{lat:.3f}°N', transform=axp.transAxes, ha='right', va='top', fontsize=7, color=colors[k])
        if k == 0:
            axp.legend(fontsize=7, ncol=3, loc='lower right')
        if k < NPROF - 1:
            axp.set_xticklabels([])
        else:
            axp.set_xlabel('Distance W→E (km)', fontsize=9)

    out = os.path.join(result_dir, 'Profile', 'profiles_combined.png')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=140, bbox_inches='tight')
    plt.close(fig)
    return out

if __name__ == '__main__':
    p = make_profile_figure(sys.argv[1] if len(sys.argv) > 1 else '.')
    print('saved', p)
