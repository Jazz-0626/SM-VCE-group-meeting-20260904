"""Rebuild SMVCE_DATA_A with 9 observations (adds LT1_Asc_POT_LOS).

Pipeline per field: crop raw -> region A @1.5" (bilinear) -> clip POT/MAI ->
linear detrend (deformation box excluded + robust) ; ALOS geometry constant ;
convex-hull AZI mask.tif ; data_information.
"""
import os, subprocess, tempfile
import numpy as np
import rasterio
from scipy.spatial import ConvexHull
from scipy.ndimage import binary_opening
from matplotlib.path import Path

BASE = os.path.dirname(os.path.abspath(__file__))
RAW = f'{BASE}/raw_renamed_0606'
RAWPOT = '/media/jazz/SSD1/TJZ/#Post/2026.04.13 小论文final/Data/0605_Rawgrd'
OUT = f'{BASE}/SMVCE_DATA_A'
A = (87.21, 28.26, 87.77, 29.07)          # W,S,E,N
RES = 0.000416666666667
DEF = (87.30, 87.66, 28.40, 29.00)        # deformation box (exclude from trend fit)
CLIP = {'S1_Asc_POT_LOS': 2.5, 'S1_Des_POT_LOS': 2.5, 'LT1_Asc_POT_LOS': 2.5,
        'LT1_Asc_POT_AZI': 2.0, 'LT1_Asc_MAI': 2.0}

# (disp_name, geom_flag, geom_prefix, raw_source_dir)
DISP = [
    ('S1_Asc_DInSAR', 1, 'S1_Asc', RAW),
    ('S1_Asc_POT_LOS', 1, 'S1_Asc', RAW),
    ('LT1_Asc_DInSAR', 1, 'LT1_Asc', RAW),
    ('LT1_Asc_POT_LOS', 1, 'LT1_Asc', RAWPOT),   # NEW 9th obs
    ('LT1_Asc_POT_AZI', 2, 'LT1_Asc', RAW),
    ('LT1_Asc_MAI', 2, 'LT1_Asc', RAW),
    ('S1_Des_DInSAR', 1, 'S1_Des', RAW),
    ('S1_Des_POT_LOS', 1, 'S1_Des', RAW),
    ('ALOS_Des_DInSAR', 1, 'ALOS_Des', RAW),
]
GEOM = ['S1_Asc', 'LT1_Asc', 'S1_Des']   # cropped from RAW; ALOS = constant

os.makedirs(OUT, exist_ok=True)
for f in os.listdir(OUT):
    if f.endswith('.tif') or f == 'data_information':
        os.remove(os.path.join(OUT, f))

def warp(src, dst):
    subprocess.run(['gdalwarp', '-overwrite', '-q', '-t_srs', 'EPSG:4326',
                    '-te', str(A[0]), str(A[1]), str(A[2]), str(A[3]), '-tr', str(RES), str(RES),
                    '-r', 'bilinear', '-of', 'GTiff', '-ot', 'Float32', '-dstnodata', 'nan', src, dst],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

print('[1] crop 9 disp -> A@1.5"')
for nm, flag, gp, srcdir in DISP:
    warp(f'{srcdir}/{nm}.disp.grd', f'{OUT}/{nm}.disp.tif')

print('[2] crop geometry (S1_Asc/LT1_Asc/S1_Des)')
for gp in GEOM:
    warp(f'{RAW}/{gp}.inc.tif', f'{OUT}/{gp}.inc.tif')
    warp(f'{RAW}/{gp}.azi.tif', f'{OUT}/{gp}.azi.tif')

print('[3] ALOS constant geometry')
with rasterio.open(f'{OUT}/ALOS_Des_DInSAR.disp.tif') as s:
    prof = s.profile; H, W = s.height, s.width; tr = s.transform
prof.update(dtype='float32', nodata=np.nan, count=1, compress='deflate')
for nm, v in [('ALOS_Des.inc.tif', 41.55), ('ALOS_Des.azi.tif', 190.07)]:
    with rasterio.open(f'{OUT}/{nm}', 'w', **prof) as d:
        d.write(np.full((H, W), v, np.float32), 1)

print('[4] clip POT/MAI')
for nm, th in CLIP.items():
    p = f'{OUT}/{nm}.disp.tif'
    with rasterio.open(p) as s:
        pf = s.profile; a = s.read(1).astype(np.float32)
    a[np.abs(a) > th] = np.nan
    with rasterio.open(p, 'w', **pf) as d:
        d.write(a, 1)

print('[5] linear detrend (deformation excluded + robust)')
lon = tr.c + (np.arange(W) + 0.5) * tr.a
lat = tr.f + (np.arange(H) + 0.5) * tr.e
LON, LAT = np.meshgrid(lon, lat)
xN = (LON - LON.mean()) / (LON.std() + 1e-9)
yN = (LAT - LAT.mean()) / (LAT.std() + 1e-9)
defbox = (LON >= DEF[0]) & (LON <= DEF[1]) & (LAT >= DEF[2]) & (LAT <= DEF[3])
for nm, flag, gp, srcdir in DISP:
    p = f'{OUT}/{nm}.disp.tif'
    with rasterio.open(p) as s:
        pf = s.profile; a = s.read(1).astype(np.float64)
    valid = np.isfinite(a)
    fit = valid & (~defbox)
    for _ in range(4):
        M = np.column_stack([np.ones(fit.sum()), xN[fit], yN[fit]])
        coef, *_ = np.linalg.lstsq(M, a[fit], rcond=None)
        pred = coef[0] + coef[1] * xN + coef[2] * yN
        r = a[fit] - pred[fit]; med = np.median(r); mad = np.median(np.abs(r - med)) * 1.4826 + 1e-9
        fit = valid & (~defbox) & (np.abs((a - pred)) < 2.5 * mad)
    a_dt = (a - (coef[0] + coef[1] * xN + coef[2] * yN)).astype(np.float32)
    with rasterio.open(p, 'w', **pf) as d:
        d.write(a_dt, 1)

print('[6] convex-hull AZI mask.tif')
def rdv(p):
    with rasterio.open(p) as s:
        a = s.read(1).astype(float); a[a == s.nodata] = np.nan
    return np.isfinite(a)
azi = rdv(f'{OUT}/LT1_Asc_POT_AZI.disp.tif') | rdv(f'{OUT}/LT1_Asc_MAI.disp.tif')
ys, xs = np.where(binary_opening(azi, iterations=3))
hv = np.column_stack([xs, ys])[ConvexHull(np.column_stack([xs, ys])).vertices]
gx, gy = np.meshgrid(np.arange(W), np.arange(H))
inside = Path(hv).contains_points(np.column_stack([gx.ravel(), gy.ravel()])).reshape(H, W)
mprof = dict(prof); mprof.update(nodata=0)
with rasterio.open(f'{OUT}/mask.tif', 'w', **mprof) as d:
    d.write(inside.astype('float32'), 1)

print('[7] data_information')
with open(f'{OUT}/data_information', 'w') as f:
    f.write('#Displacement_Measurements\tLOS_AZI_ENU\tInc_Angle\tAzi_Angle\tLeft_or_Right_looking\n')
    for nm, flag, gp, srcdir in DISP:
        f.write(f'{nm}.disp.tif\t{flag}\t{gp}.inc.tif\t{gp}.azi.tif\t1\n')

# consistency + summary
print('\n[check] grid + valid%')
ref = None
for nm, flag, gp, srcdir in DISP:
    with rasterio.open(f'{OUT}/{nm}.disp.tif') as s:
        a = s.read(1).astype(float); a[a == s.nodata] = np.nan
        k = (s.width, s.height)
        ref = ref or k
        flagc = 'LOS' if flag == 1 else 'AZI'
        print(f'  {nm:18s}[{flagc}] {s.width}x{s.height} valid={100*np.isfinite(a).mean():4.1f}% range=[{np.nanmin(a):+.2f},{np.nanmax(a):+.2f}]')
print(f'  mask solve pixels = {int(inside.sum())} ({100*inside.mean():.1f}%)')
print(f'  grid {ref[0]}x{ref[1]}  -> 9 obs ready in {OUT}')
