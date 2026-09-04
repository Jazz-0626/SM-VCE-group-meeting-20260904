"""Re-plot 3D_SMVCE.png from an existing result folder (no re-solving).

Improvements:
  - per-component symmetric colorbar limits (so the vertical >2 m is not clipped)
  - edge cleanup for N-S / U (erode boundary ring + physical outlier clip)
  - denser, longer, higher-contrast (black) horizontal quiver, kept off the mask edge

Usage: python replot_3d.py <result_dir>
"""
import os, sys
os.environ.pop('PROJ_DATA', None); os.environ.pop('PROJ_LIB', None); os.environ.pop('GDAL_DATA', None)
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import rasterio
from scipy.ndimage import distance_transform_edt, binary_erosion

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from smvce_code.plotting import apply_publication_style, getfig, coortick, plot_fault, clbtitle

res = sys.argv[1] if len(sys.argv) > 1 else '.'
g = os.path.join(res, 'geotiff')

def rd(name):
    with rasterio.open(os.path.join(g, name)) as s:
        a = s.read(1).astype(float)
        if s.nodata is not None:
            a[a == s.nodata] = np.nan
        tr = s.transform
    return a, tr

e, tr = rd('enu_east.tif'); n, _ = rd('enu_north.tif'); u, _ = rd('enu_up.tif')
H, W = e.shape
coor = {'corner_lon': tr.c, 'corner_lat': tr.f, 'post_lon': tr.a, 'post_lat': tr.e,
        'width': W, 'nlines': H}
fault = 0

# ---- edge cleanup: erode the solved region + physical outlier clip ----
solved = np.isfinite(e) & np.isfinite(n) & np.isfinite(u)
core = binary_erosion(solved, iterations=15)         # drop ~15 px boundary ring (aggressive)
PHYS = {'e': 1.5, 'n': 1.0, 'u': 3.0}                 # physical caps (m); beyond = artifact
def clean(a, cap):
    b = np.where(core, a, np.nan)
    b[np.abs(b) > cap] = np.nan
    return b
e = clean(e, PHYS['e']); n = clean(n, PHYS['n']); u = clean(u, PHYS['u'])
enu = np.stack([e, n, u], axis=2)

# ---- per-component symmetric color limits (cover the real range, e.g. U ~ 2.4) ----
# per-component fixed colour limits: E ±0.8, N ±0.8, U ±2 (edit here)
dl = np.array([[-0.8, 0.8], [-0.8, 0.8], [-2.0, 2.0]])
print(f"colorbar limits: E/N ±0.8, U ±2.0 m  (U real min={np.nanmin(u):.2f})")

apply_publication_style()
lgd = ['(a) E-W', '(b) N-S', '(c) Vertical']
fig, axs, cbs = getfig(enu, dl, flag_clb=True, lgdstr=lgd)
for idx, ax in enumerate(axs):
    coortick(ax, coor, dxy=(0.2, 0.2))
    if idx < len(cbs):
        clbtitle(cbs[idx], '[m]')
    plot_fault(ax, fault, coor)

# ---- denser, longer, high-contrast quiver on the vertical panel ----
nv = 16
ri = np.round(np.linspace(0, H - 1, nv)).astype(int)
ci = np.round(np.linspace(0, W - 1, nv)).astype(int)
CC, RR = np.meshgrid(ci, ri)
valid2d = np.isfinite(enu[:, :, 0]) & np.isfinite(enu[:, :, 1])
edge = distance_transform_edt(valid2d)
margin = 0.6 * (H / nv)
keep = valid2d[RR, CC] & (edge[RR, CC] > margin)
E = enu[RR, CC, 0]; N = enu[RR, CC, 1]
keep &= np.isfinite(E) & np.isfinite(N)
if np.any(keep):
    ax_v = axs[2]
    QSCALE = 3.0   # smaller -> longer arrows (more visible length variation)
    ax_v.quiver(CC[keep], RR[keep], E[keep], N[keep],
                color='black', alpha=0.9, scale=QSCALE, scale_units='width',
                width=0.005, headwidth=4.5, headlength=5.5, headaxislength=4.5, zorder=6)
    rx, ry = 0.10 * W, 0.93 * H
    ax_v.quiver(rx, ry, 0.5, 0, color='black', alpha=1.0, scale=QSCALE, scale_units='width',
                width=0.005, headwidth=4.5, headlength=5.5, headaxislength=4.5, zorder=7)
    ax_v.text(rx, ry + 0.035 * H, '0.5 m', fontsize=9, color='black', ha='left', va='bottom',
              fontweight='bold', bbox=dict(boxstyle='round,pad=0.15', facecolor='white', edgecolor='none', alpha=0.8))

out = os.path.join(res, '3D_SMVCE.png')
fig.savefig(out, dpi=150, bbox_inches='tight')
plt.close(fig)
print(f"saved {out}  (quiver nv={nv}, kept {int(keep.sum())} arrows, scale={QSCALE})")
