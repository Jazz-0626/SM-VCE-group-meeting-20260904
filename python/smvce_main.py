"""
SM-VCE Main Script - Python Version
Calculate 3-D surface displacements of the 8 January 2022 Mw6.7 Menyuan earthquake.

Based on the Strain Model and Variance Component Estimation (SM-VCE) method.
Original MATLAB code by Dr. Jihong Liu, liujihong@csu.edu.cn
Converted to Python.

Usage:
    python smvce_main.py
"""
import os
import sys
import time
import shutil
import json
import argparse

# Windows/Conda compatibility: this environment contains both expat.dll and
# libexpat.dll.  Preload the environment-local DLL before Rasterio indirectly
# imports pyexpat, otherwise Windows may bind pyexpat to an incompatible copy
# earlier on PATH.  This only fixes DLL resolution; it does not alter SM-VCE.
if sys.platform == 'win32':
    try:
        import ctypes
        _conda_dll_dir = os.path.join(sys.prefix, 'Library', 'bin')
        _conda_dll_dir_handle = os.add_dll_directory(_conda_dll_dir)
        _expat_dll = os.path.join(_conda_dll_dir, 'expat.dll')
        if os.path.isfile(_expat_dll):
            ctypes.WinDLL(_expat_dll)
            import pyexpat  # noqa: F401 - bind before Rasterio imports XML support
    except OSError:
        pass

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for saving figures
import matplotlib.pyplot as plt

# Add code directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# --- PROJ/GDAL data isolation -------------------------------------------------
# rasterio is installed as a manylinux pip wheel that bundles its own
# libproj (PROJ 9.4.1, db layout >= 1.3) and rasterio/proj_data/proj.db.
# When the conda env is `activate`d, proj4-activate.sh exports
# PROJ_DATA -> $CONDA_PREFIX/share/proj (conda PROJ 9.3.1, db layout 1.2),
# which the bundled libproj then rejects with:
#   "proj.db ... LAYOUT.VERSION.MINOR = 2 whereas a number >= 3 is expected".
# Clearing these env vars before rasterio is imported makes rasterio fall back
# to its own bundled proj_data, while any subprocess GMT/GDAL still use their
# system defaults. (System/conda/pip ship three different PROJ versions, so a
# single global PROJ_DATA cannot satisfy all of them.)
for _v in ('PROJ_DATA', 'PROJ_LIB', 'GDAL_DATA'):
    os.environ.pop(_v, None)
# -----------------------------------------------------------------------------

from smvce_code.smvce_readdispdata import smvce_readdispdata
from smvce_code.wls3d import wls3d
from smvce_code.smvce_solve3d import smvce_solve3d
from smvce_code.export_results import export_smvce_results
from smvce_code.plotting import (
    apply_publication_style, getfig, coortick, plot_fault, plot2d_vector,
    clbtitle, save_figure, compute_dls, select_sparse_vectors
)
from smvce_code.utils import sub2lonlat


# =====================================================
# Figure / output configuration (edit here)
# =====================================================
# Per-component figure colour limits, metres: [[E_lo,E_hi],[N_lo,N_hi],[U_lo,U_hi]].
DL_ENU = [[-0.8, 0.8], [-0.8, 0.8], [-2.0, 2.0]]
# Physical cap on the solved N-S amplitude (m). Pixels whose |N-S| exceeds this are
# set to NaN, removing ill-conditioned spikes at pixels lacking azimuth coverage.
# Set to None to disable.
NS_OUTPUT_CLIP = 1.0


def main():
    parser = argparse.ArgumentParser(description='SM-VCE 3D Deformation Solver')
    parser.add_argument('--windowsize', '-w', type=int, default=50,
                        help='SM-VCE strain-model window size (default: 50)')
    parser.add_argument('--tag', type=str, default='',
                        help='Optional tag appended to result folder name')
    parser.add_argument('--data-dir', '-d', type=str, default=None,
                        help='Path to the SM-VCE input data folder (must contain a '
                             'data_information file). Accepts an absolute path, or a '
                             'name/relative path resolved against the script directory '
                             'and the current working directory. '
                             'Default: auto-detect a folder named "SMVCE_DATA".')
    parser.add_argument('--no-wls', action='store_true',
                        help='Skip the classical WLS inversion and run SM-VCE only '
                             '(faster).')
    parser.add_argument('--exclude-obs', default='',
                        help='Comma-separated 1-based observation indices to exclude, e.g. 3,6.')
    parser.add_argument('--output-dir', default='.',
                        help='Parent directory for the timestamped result folder.')
    parser.add_argument('--no-cuda', action='store_true', help='Disable the CUDA fast path.')
    parser.add_argument('--preliminary-gpu-no-fault-separation', action='store_true',
                        help='Generate a full-scene GPU preview without cross-fault '
                             'window separation. This is a preview-only acceleration; '
                             'the fault trace is still plotted.')
    parser.add_argument('--no-grd', action='store_true',
                        help='Skip GMT GRD export; GeoTIFF output is still written.')
    parser.add_argument('--ns-output-clip', type=float, default=None,
                        help='Case-specific N-S absolute clip in metres; disabled by default.')
    parser.add_argument('--figure-erode', type=int, default=0,
                        help='Case-specific figure-only edge erosion in pixels; default 0.')
    parser.add_argument('--cuda-rows-per-batch', type=int, default=None,
                        help='Override automatic CUDA row batch size; useful on 8 GB GPUs.')
    parser.add_argument('--roi', nargs=4, type=int, metavar=('R1', 'R2', 'C1', 'C2'),
                        help='Optional 1-based inclusive crop: row_start row_end col_start col_end.')
    args = parser.parse_args()

    # =====================================================
    # SM-VCE Configuration
    # =====================================================

    # Data path - point to the original MATLAB data folder
    script_dir = os.path.dirname(os.path.abspath(__file__))
    matlab_dir = os.path.join(os.path.dirname(script_dir),
                              'SMVCE_DEMO_20250814_Menyuan Earthquake')

    # Resolve the input data folder.
    # Priority: --data-dir (explicit) > local SMVCE_DATA > MATLAB demo SMVCE_DATA.
    if args.data_dir:
        candidates = [
            args.data_dir,
            os.path.join(script_dir, args.data_dir),
            os.path.join(os.getcwd(), args.data_dir),
        ]
        pdata = next((c for c in candidates if os.path.isdir(c)), None)
        if pdata is None:
            print(f"错误：指定的数据目录不存在：{args.data_dir}")
            print("已尝试以下位置：")
            for c in candidates:
                print(f"  - {c}")
            sys.exit(1)
    elif os.path.exists(os.path.join(script_dir, 'SMVCE_DATA')):
        pdata = os.path.join(script_dir, 'SMVCE_DATA')
    elif os.path.exists(os.path.join(matlab_dir, 'SMVCE_DATA')):
        pdata = os.path.join(matlab_dir, 'SMVCE_DATA')
    else:
        print("错误：未找到数据目录！")
        print("请用 --data-dir/-d 显式指定，或将数据放入默认的 SMVCE_DATA 目录。")
        print(f"预期默认位置：{os.path.join(script_dir, 'SMVCE_DATA')}")
        sys.exit(1)

    pdata = os.path.abspath(pdata)
    # Sanity check: the chosen folder must contain a data_information file.
    if not os.path.isfile(os.path.join(pdata, 'data_information')):
        print(f"错误：数据目录缺少 data_information 文件：{pdata}")
        sys.exit(1)

    print(f"数据来源：{pdata}")
    apply_publication_style()

    # Read displacement data
    print("正在读取形变观测数据...", flush=True)
    data, inc, azi, losazienu, leftorright, coor, dem, mask, fault, datainfo = \
        smvce_readdispdata(pdata)

    if args.roi:
        r1, r2, c1, c2 = args.roi
        if not (1 <= r1 <= r2 <= data.shape[0] and 1 <= c1 <= c2 <= data.shape[1]):
            raise SystemExit(f'--roi outside raster bounds 1..{data.shape[0]}, 1..{data.shape[1]}')
        rs, cs = slice(r1 - 1, r2), slice(c1 - 1, c2)
        data = data[rs, cs, :]
        inc = inc[rs, cs, :]
        azi = azi[rs, cs, :]
        dem = dem[rs, cs]
        mask = mask[rs, cs]
        coor['corner_lon'] += (c1 - 1) * coor['post_lon']
        coor['corner_lat'] += (r1 - 1) * coor['post_lat']
        coor['nlines'], coor['width'] = data.shape[:2]
        print(f'ROI：rows {r1}:{r2}, cols {c1}:{c2} -> {data.shape[0]}x{data.shape[1]}')

    # =====================================================
    # Parameters
    # =====================================================
    # Window size for the strain model (CLI: -w / --windowsize, default 50).
    windowsize = args.windowsize
    print(f"窗口大小 windowsize = {windowsize}")
    # 0: have azimuth obs -> 3D; 1: only LOS -> 2D (E-W and Vertical)
    flag_if_2D = 0

    # Strain model dimension (2 recommended for real cases)
    fsmpara = 2

    # Inter-weight flag (relative weighting within window)
    flag_interWeight = 0

    # Adaptive window size flag
    flag_adpws = 0

    # SMAD algorithm flag (adaptive fault removal)
    flag_smad = 0

    # CUDA acceleration flag for the fixed-window default SM-VCE path
    use_cuda = not args.no_cuda

    # Number of image rows processed per CUDA batch (auto-adapt for large windows)
    # Memory per batch ~ rows_batch * cols * ws^2 * data_num * 8 bytes
    cuda_rows_per_batch = (args.cuda_rows_per_batch if args.cuda_rows_per_batch is not None
                           else max(2, int(16 * (40 / windowsize) ** 2)))  # ws40=16, ws60=7, ws80=4
    if cuda_rows_per_batch < 1:
        raise SystemExit('--cuda-rows-per-batch must be >= 1')

    # Whether to run the classical WLS inversion (skip with --no-wls for speed).
    run_wls = not args.no_wls

    # Whether to run the SM-VCE inversion.
    run_smvce = True

    # Indices of observations NOT to use (0-based)
    try:
        indnotuse = [int(x.strip()) - 1 for x in args.exclude_obs.split(',') if x.strip()]
    except ValueError as exc:
        raise SystemExit(f'--exclude-obs must contain comma-separated integers: {exc}')
    if any(i < 0 or i >= data.shape[2] for i in indnotuse):
        raise SystemExit(f'--exclude-obs outside valid 1..{data.shape[2]} range')

    # Fixed ENU colorbar limits for 3D_WLS.png and 3D_SMVCE.png.
    dl_enu_plot = np.array([
        [-1.5, 1.5],
        [-1.5, 1.5],
        [-1.5, 1.5],
    ], dtype=float)

    vector_max_count_plot = 36

    # =====================================================
    # Show input data
    # =====================================================
    datanms = [info[0] for info in datainfo]
    data = data.astype(float)
    data_plot = data.copy()
    data_plot[data_plot == 0] = np.nan

    print("正在绘制输入观测数据图 DATA.png ...", flush=True)
    # Per-observation color limits: LOS=[-2,2], AZI=[-1,1]
    dl_data = np.zeros((data_plot.shape[2], 2))
    for idx in range(data_plot.shape[2]):
        if losazienu[idx] == 1:  # LOS
            dl_data[idx] = [-2.0, 2.0]
        else:                     # AZI
            dl_data[idx] = [-1.0, 1.0]
    fig, axs, cbs = getfig(data_plot, dl_data, flag_clb=True, lgdstr=datanms)
    for idx, ax in enumerate(axs):
        coortick(ax, coor, dxy=(0.2, 0.2))
        if idx < len(cbs):
            clbtitle(cbs[idx], '[m]')
        plot_fault(ax, fault, coor)
    save_figure(fig, 'DATA.png')
    plt.close(fig)

    # =====================================================
    # Remove unused observations
    # =====================================================
    if len(indnotuse) > 0:
        data = np.delete(data, indnotuse, axis=2)
        inc = np.delete(inc, indnotuse, axis=2)
        azi = np.delete(azi, indnotuse, axis=2)
        losazienu = np.delete(losazienu, indnotuse)
        leftorright = np.delete(leftorright, indnotuse)
        for idx in sorted(indnotuse, reverse=True):
            del datainfo[idx]
    datanms = [info[0] for info in datainfo]

    # Replace NaN with 0 for processing (0 is excluded by solver as "no data")
    data[np.isnan(data)] = 0

    # Prepare data dictionary (equivalent to MATLAB's save/load)
    solver_fault = 0 if args.preliminary_gpu_no_fault_separation else fault
    if args.preliminary_gpu_no_fault_separation:
        print('整景预览模式：GPU 求解不执行断层跨侧窗口剔除；断层仍叠加在结果图中。')

    data_dict = {
        'data': data, 'inc': inc, 'azi': azi,
        'losazienu': losazienu, 'flag_if_2D': flag_if_2D,
        'leftorright': leftorright, 'coor': coor,
        'dem': dem, 'mask': mask, 'fault': solver_fault,
        'windowsize': windowsize, 'fsmpara': fsmpara,
        'flag_interWeight': flag_interWeight,
        'flag_adpws': flag_adpws, 'flag_smad': flag_smad,
        'use_cuda': use_cuda,
        'cuda_rows_per_batch': cuda_rows_per_batch
    }

    lgdstr_3d = ['(a) E-W', '(b) N-S', '(c) Vertical']
    Result_wls = None
    Result_smvce = None

    # =====================================================
    # WLS: Classical Weighted Least Squares
    # =====================================================
    if run_wls:
        Result_wls = wls3d(data_dict)

        fig, axs, cbs = getfig(Result_wls['enu'], dl_enu_plot, flag_clb=True, lgdstr=lgdstr_3d)
        for idx, ax in enumerate(axs):
            coortick(ax, coor, dxy=(0.2, 0.2))
            if idx < len(cbs):
                clbtitle(cbs[idx], '[m]')
            plot_fault(ax, fault, coor)

        row, col, _ = Result_wls['enu'].shape
        xyen = select_sparse_vectors(Result_wls['enu'], max_count=vector_max_count_plot)
        plot2d_vector(axs[2], xyen, field_shape=(row, col))

        save_figure(fig, '3D_WLS.png')
        plt.close(fig)
    else:
        print('跳过 WLS 反演。')

    # =====================================================
    # SM-VCE: Strain Model + Variance Component Estimation
    # =====================================================
    if run_smvce:
        Result_smvce = smvce_solve3d(data_dict)

        # Physical output clip on N-S: remove ill-conditioned spikes (pixels lacking
        # azimuth coverage) so the exported field and figures are clean.
        if args.ns_output_clip is not None:
            _ns = Result_smvce['enu'][:, :, 1]
            _nbad = int(np.sum(np.abs(_ns) > args.ns_output_clip))
            _ns[np.abs(_ns) > args.ns_output_clip] = np.nan
            Result_smvce['enu'][:, :, 1] = _ns
            print(f"N-S 震例阈值 |N-S|>{args.ns_output_clip}m -> NaN：剔除 {_nbad} 个像元")

        # Save raster results
        timestamp = time.strftime('%Y%m%d%H%M%S')
        tag = f'_ws{windowsize}'
        if args.tag:
            tag += f'_{args.tag}'
        output_root = os.path.join(os.path.abspath(args.output_dir), f'Result_SMVCE_{timestamp}{tag}')
        os.makedirs(output_root, exist_ok=True)
        write_grd = (not args.no_grd) and shutil.which('gmt') is not None
        if not write_grd:
            print('未检测到 GMT 或已指定 --no-grd：跳过 GRD，仅导出 GeoTIFF。')
        export_info = export_smvce_results(Result_smvce, coor, datainfo, output_root, write_grd=write_grd)
        print(f"GeoTIFF 结果已导出至：{export_info['geotiff_dir']}")
        if write_grd:
            print(f"GRD 结果已导出至：{export_info['grd_dir']}")
        print(f"结果说明文件 README 已写入：{export_info['readme_path']}")

        # =====================================================
        # Plot SM-VCE 3D Displacement Results
        # =====================================================
        # --- display cleanup: aggressive edge removal so N-S/U boundary artifacts
        #     do not show (erode the solved region + physical outlier clip).
        #     Only the FIGURE is cleaned; the exported GeoTIFF/GRD keep the raw result.
        from scipy.ndimage import distance_transform_edt, binary_erosion
        _enu_raw = Result_smvce['enu']
        _row, _col, _ = _enu_raw.shape
        _solved = (np.isfinite(_enu_raw[:, :, 0]) & np.isfinite(_enu_raw[:, :, 1])
                   & np.isfinite(_enu_raw[:, :, 2]))
        _core = binary_erosion(_solved, iterations=args.figure_erode) if args.figure_erode > 0 else _solved
        _phys = [1.5, 1.0, 3.0]  # E, N, U physical caps (m); beyond -> artifact
        _enu = _enu_raw.copy()
        for _k in range(3):
            _lay = _enu[:, :, _k]
            _lay[~_core] = np.nan
            _lay[np.abs(_lay) > _phys[_k]] = np.nan
            _enu[:, :, _k] = _lay

        # per-component figure colour limits (configured at top: DL_ENU)
        _dl = np.array(DL_ENU)

        with open(os.path.join(output_root, 'run_config.json'), 'w', encoding='utf-8') as _cf:
            json.dump({
                'data_dir': pdata,
                'windowsize': windowsize,
                'excluded_observations_1based': [i + 1 for i in indnotuse],
                'fsmpara': fsmpara,
                'use_cuda_requested': use_cuda,
                'cuda_rows_per_batch': cuda_rows_per_batch,
                'preliminary_gpu_no_fault_separation': args.preliminary_gpu_no_fault_separation,
                'ns_output_clip_m': args.ns_output_clip,
                'figure_erode_pixels': args.figure_erode,
                'write_grd': write_grd,
                'run_wls': run_wls,
                'roi_1based_inclusive': args.roi,
            }, _cf, ensure_ascii=False, indent=2)

        fig, axs, cbs = getfig(_enu, _dl, flag_clb=True, lgdstr=lgdstr_3d)
        for idx, ax in enumerate(axs):
            coortick(ax, coor, dxy=(0.2, 0.2))
            if idx < len(cbs):
                clbtitle(cbs[idx], '[m]')
            plot_fault(ax, fault, coor)

        # Horizontal vectors on the vertical subplot: denser, longer, high-contrast
        # (black), kept off the mask edge.
        _nv = 16
        _ri_samp = np.round(np.linspace(0, _row - 1, _nv)).astype(int)
        _ci_samp = np.round(np.linspace(0, _col - 1, _nv)).astype(int)
        _CC, _RR = np.meshgrid(_ci_samp, _ri_samp)
        _E = _enu[_RR, _CC, 0]
        _N = _enu[_RR, _CC, 1]
        _valid2d = np.isfinite(_enu[:, :, 0]) & np.isfinite(_enu[:, :, 1])
        _edge_dist = distance_transform_edt(_valid2d)
        _margin_px = 0.6 * (_row / _nv)
        _valid_q = (_valid2d[_RR, _CC] & (_edge_dist[_RR, _CC] > _margin_px)
                    & np.isfinite(_E) & np.isfinite(_N))
        if np.any(_valid_q):
            _ax_v = axs[2]
            _qscale = 3.0  # smaller -> longer arrows (clearer length variation)
            _ax_v.quiver(
                _CC[_valid_q], _RR[_valid_q],
                _E[_valid_q], _N[_valid_q],
                color='black', alpha=0.9, scale=_qscale, scale_units='width',
                width=0.005, headwidth=4.5, headlength=5.5, headaxislength=4.5,
                zorder=6
            )
            _ref_mag = 0.5
            _rx, _ry = 0.10 * _col, 0.93 * _row
            _ax_v.quiver(
                _rx, _ry, _ref_mag, 0,
                color='black', alpha=1.0, scale=_qscale, scale_units='width',
                width=0.005, headwidth=4.5, headlength=5.5, headaxislength=4.5,
                zorder=7
            )
            _ax_v.text(_rx, _ry + 0.035 * _row, f'{_ref_mag} m', fontsize=9,
                      color='black', ha='left', va='bottom', fontweight='bold',
                      bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                               edgecolor='none', alpha=0.8))

        save_figure(fig, '3D_SMVCE.png')
        plt.close(fig)

        # =====================================================
        # Multi-profile analysis (combined figure: E-W map + 10 ENU profiles)
        # =====================================================
        print("正在绘制近场多剖面组合图 ...", flush=True)
        try:
            from plot_profiles import make_profile_figure
            _pf = make_profile_figure(output_root)
            print(f"  已保存组合剖面图 -> {_pf}")
        except Exception as _pe:
            print(f"  [WARN] 剖面图生成失败: {_pe}")

        # =====================================================
        # Compute accuracy statistics in far-field area
        # =====================================================
        row, col, _ = Result_smvce['enu'].shape
        ii_far = np.arange(int(round(0.8 * row)), row)
        jj_far = np.arange(int(round(0.8 * col)), col)
        enustd = Result_smvce['enu'][np.ix_(ii_far, jj_far, [0, 1, 2])]
        enustd_flat = enustd.reshape(-1, 3)
        enustdv = np.nanstd(enustd_flat, axis=0)
        print(f'\n远场精度（std）：E-W = {enustdv[0]:.4f} m，N-S = {enustdv[1]:.4f} m，Vertical = {enustdv[2]:.4f} m')

        recsub = np.array([
            [ii_far[0], jj_far[0]],
            [ii_far[0], jj_far[-1]],
            [ii_far[-1], jj_far[-1]],
            [ii_far[-1], jj_far[0]],
            [ii_far[0], jj_far[0]]
        ])
        recll = sub2lonlat(recsub, coor)
        print(f'远场统计区域：lon = [{recll[0, 0]:.3f}, {recll[2, 0]:.3f}]，'
              f'lat = [{recll[2, 1]:.3f}, {recll[0, 1]:.3f}]')

        # =====================================================
        # Plot Standard Deviation
        # =====================================================
        log10_std_enu = np.log10(np.maximum(Result_smvce['var']['enu'], 1e-15))
        dl_std = compute_dls(log10_std_enu, symmetric=False)
        lgdstr_std = ['(a) E-W Std', '(b) N-S Std', '(c) Vertical Std']
        fig, axs, cbs = getfig(log10_std_enu, dl_std, flag_clb=True, lgdstr=lgdstr_std)
        for idx, ax in enumerate(axs):
            coortick(ax, coor, dxy=(0.2, 0.2))
            if idx < len(cbs):
                clbtitle(cbs[idx], r'$\log_{10}(\mathrm{std}\,[m])$')
            plot_fault(ax, fault, coor)

        save_figure(fig, '3D_std_SMVCE.png')
        plt.close(fig)

        # =====================================================
        # Copy all result figures into the Result folder
        # =====================================================
        _fig_names = ['DATA.png', '3D_SMVCE.png', '3D_std_SMVCE.png']
        if run_wls:
            _fig_names.append('3D_WLS.png')
        print(f"\n正在将结果图复制到 {output_root}/")
        for _fn in _fig_names:
            if os.path.exists(_fn):
                shutil.copy2(_fn, os.path.join(output_root, _fn))
                print(f"  -> {output_root}/{_fn}")
    else:
        print('跳过 SM-VCE 反演。')

    if not run_wls and not run_smvce:
        print('未选择任何反演方法，仅生成了输入观测数据图。')

    print("\n" + "=" * 60)
    print("全部处理完成！")
    if run_smvce:
        print(f"结果保存于：{output_root}/")
    print("=" * 60)


if __name__ == '__main__':
    main()
