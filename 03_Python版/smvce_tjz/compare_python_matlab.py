"""
compare_python_matlab.py
========================

逐像元对比 Python 版与 MATLAB 版 SM-VCE 输出的三维形变场（E-W / N-S / Vertical）。

设计的指标：
    1. 有效像元覆盖率（两边都有效的像元 / 总像元）
    2. 均值差 (Mean diff = Python − MATLAB)
    3. 标准差 (Std of diff)
    4. RMSE √(mean((py − mat)²))
    5. MAE  mean(|py − mat|)
    6. Max |diff|
    7. Pearson 相关系数 r
    8. 落在 ±0.05 m / ±0.10 m / ±0.20 m 内的像元占比

可视化：
    1. Comparison_<component>.png  3 子图（Python | MATLAB | Diff）+ colorbar
    2. Scatter_<component>.png     散点 Python vs MATLAB（含拟合线、1:1 线）
    3. Histogram_diff_<component>.png  差值分布直方图
    4. Summary.png                  总览（E/N/U × 三类图）
    5. metrics.txt                  数值指标全文文本
    6. metrics.csv                  数值指标 CSV，便于后续脚本读取

用法：
    conda run -n smvce_tiff python compare_python_matlab.py
    conda run -n smvce_tiff python compare_python_matlab.py \\
        --python-dir Result_SMVCE_20260528211710_ws55/geotiff \\
        --matlab-dir "/home/.../Result_SMVCE_20260528233509_matlab_geotiff"

缺省自动定位：
    --python-dir 选 Result_SMVCE_*/geotiff/ 中时间戳最新的一个
    --matlab-dir 选 /home/jazz/tools/SMVCE_DEMO_20250814_Menyuan Earthquake/
                  下 Result_SMVCE_*_matlab_geotiff/ 中时间戳最新的一个
"""
import argparse
import csv
import glob
import os
import sys
from datetime import datetime

import numpy as np
import rasterio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

# 配置中文字体（系统装有 Noto Sans CJK JP），避免标题里的中文字符变成方框
plt.rcParams['font.sans-serif'] = ['Noto Sans CJK JP', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False


COMPONENTS = [
    ('east', '东向 E-W'),
    ('north', '北向 N-S'),
    ('up', '垂向 Vertical'),
]


# -----------------------------------------------------------------------------
# 输入定位
# -----------------------------------------------------------------------------
def auto_locate_python_dir() -> str:
    here = os.path.dirname(os.path.abspath(__file__))
    candidates = sorted(glob.glob(os.path.join(here, 'Result_SMVCE_*/geotiff')))
    if not candidates:
        print('错误：未在 Python 项目目录下找到 Result_SMVCE_*/geotiff/')
        sys.exit(1)
    return candidates[-1]


def auto_locate_matlab_dir() -> str:
    matlab_root = '/home/jazz/tools/SMVCE_DEMO_20250814_Menyuan Earthquake'
    candidates = sorted(glob.glob(os.path.join(matlab_root, 'Result_SMVCE_*_matlab_geotiff')))
    if not candidates:
        print(f'错误：未找到 {matlab_root}/Result_SMVCE_*_matlab_geotiff/')
        print('请先在 MATLAB 目录下运行 mat_to_geotiff.py 提取脚本。')
        sys.exit(1)
    return candidates[-1]


def load_layer(directory: str, name: str) -> tuple:
    """读取单个 GeoTIFF 图层，返回 (data, transform, shape)。"""
    path = os.path.join(directory, f'enu_{name}.tif')
    if not os.path.exists(path):
        raise FileNotFoundError(f'未找到 {path}')
    with rasterio.open(path) as s:
        arr = s.read(1).astype(np.float64)
        return arr, s.transform, s.shape


# -----------------------------------------------------------------------------
# 指标计算
# -----------------------------------------------------------------------------
def compute_metrics(py: np.ndarray, mat: np.ndarray) -> dict:
    """对单个分量计算所有指标。"""
    py_valid = np.isfinite(py) & (py != 0)
    mat_valid = np.isfinite(mat) & (mat != 0)
    both = py_valid & mat_valid
    diff = py - mat
    diff_valid = diff[both]

    metrics = {
        'n_total': int(py.size),
        'n_python_valid': int(py_valid.sum()),
        'n_matlab_valid': int(mat_valid.sum()),
        'n_both_valid': int(both.sum()),
        'coverage_pct': float(100.0 * both.sum() / py.size),
    }

    if both.sum() == 0:
        for k in ('mean_diff', 'std_diff', 'rmse', 'mae',
                  'max_abs_diff', 'corr', 'pct_within_5cm',
                  'pct_within_10cm', 'pct_within_20cm',
                  'py_mean', 'mat_mean', 'py_std', 'mat_std'):
            metrics[k] = float('nan')
        return metrics

    metrics['py_mean'] = float(np.mean(py[both]))
    metrics['mat_mean'] = float(np.mean(mat[both]))
    metrics['py_std'] = float(np.std(py[both]))
    metrics['mat_std'] = float(np.std(mat[both]))
    metrics['mean_diff'] = float(np.mean(diff_valid))
    metrics['std_diff'] = float(np.std(diff_valid))
    metrics['rmse'] = float(np.sqrt(np.mean(diff_valid ** 2)))
    metrics['mae'] = float(np.mean(np.abs(diff_valid)))
    metrics['max_abs_diff'] = float(np.max(np.abs(diff_valid)))

    # Pearson 相关
    if py[both].std() > 0 and mat[both].std() > 0:
        metrics['corr'] = float(np.corrcoef(py[both], mat[both])[0, 1])
    else:
        metrics['corr'] = float('nan')

    metrics['pct_within_5cm'] = float(100.0 * (np.abs(diff_valid) < 0.05).sum() / diff_valid.size)
    metrics['pct_within_10cm'] = float(100.0 * (np.abs(diff_valid) < 0.10).sum() / diff_valid.size)
    metrics['pct_within_20cm'] = float(100.0 * (np.abs(diff_valid) < 0.20).sum() / diff_valid.size)
    return metrics


# -----------------------------------------------------------------------------
# 可视化
# -----------------------------------------------------------------------------
def plot_comparison(py: np.ndarray, mat: np.ndarray, transform,
                    component_key: str, component_label: str,
                    output_path: str, vmin=-1.5, vmax=1.5, diff_lim=0.5):
    """3 子图：Python | MATLAB | Diff（同色标）。"""
    row, col = py.shape
    extent = [
        transform.c,                                 # left lon
        transform.c + transform.a * col,             # right lon
        transform.f + transform.e * row,             # bottom lat
        transform.f                                  # top lat
    ]

    fig, axes = plt.subplots(1, 3, figsize=(18, 7))
    norm_data = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
    norm_diff = TwoSlopeNorm(vmin=-diff_lim, vcenter=0, vmax=diff_lim)

    panels = [
        (axes[0], py, f'(a) Python {component_label}', norm_data, 'RdBu_r'),
        (axes[1], mat, f'(b) MATLAB {component_label}', norm_data, 'RdBu_r'),
        (axes[2], py - mat, f'(c) Diff (Python − MATLAB)', norm_diff, 'RdBu_r'),
    ]
    for ax, d, title, norm, cmap in panels:
        d_show = np.where(np.isfinite(d) & (d != 0), d, np.nan)
        im = ax.imshow(d_show, extent=extent, cmap=cmap, norm=norm,
                       aspect='auto', origin='upper', interpolation='nearest')
        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.set_xlabel('Longitude (°E)')
        if ax is axes[0]:
            ax.set_ylabel('Latitude (°N)')
        plt.colorbar(im, ax=ax, shrink=0.8, label='[m]')

    fig.suptitle(f'Python vs MATLAB — {component_label}', fontsize=14, fontweight='bold')
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  已保存：{output_path}')


def plot_scatter(py: np.ndarray, mat: np.ndarray,
                 component_label: str, output_path: str,
                 max_points: int = 50000):
    """散点图：MATLAB 横轴 vs Python 纵轴，含 1:1 线与最小二乘拟合。"""
    both = np.isfinite(py) & np.isfinite(mat) & (py != 0) & (mat != 0)
    py_v = py[both]
    mat_v = mat[both]
    if py_v.size == 0:
        return
    # 大数据下采样
    if py_v.size > max_points:
        idx = np.random.default_rng(0).choice(py_v.size, size=max_points, replace=False)
        py_v = py_v[idx]
        mat_v = mat_v[idx]

    lo = min(mat_v.min(), py_v.min())
    hi = max(mat_v.max(), py_v.max())
    pad = 0.05 * (hi - lo)
    lo -= pad
    hi += pad

    # 最小二乘 y = a + b·x
    A = np.vstack([np.ones_like(mat_v), mat_v]).T
    coef, *_ = np.linalg.lstsq(A, py_v, rcond=None)
    a, b = float(coef[0]), float(coef[1])
    r = float(np.corrcoef(mat_v, py_v)[0, 1])

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.hexbin(mat_v, py_v, gridsize=80, cmap='viridis', mincnt=1, bins='log')
    ax.plot([lo, hi], [lo, hi], 'k--', lw=1.5, label='1:1 line', alpha=0.7)
    ax.plot([lo, hi], [a + b * lo, a + b * hi], 'r-', lw=1.2,
            label=f'Fit: y = {b:.4f}·x + {a:.4f}  (r={r:.4f})')
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel('MATLAB [m]')
    ax.set_ylabel('Python [m]')
    ax.set_title(f'Python vs MATLAB scatter — {component_label}', fontweight='bold')
    ax.legend(loc='upper left')
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  已保存：{output_path}')


def plot_diff_histogram(py: np.ndarray, mat: np.ndarray,
                        component_label: str, output_path: str):
    """差值直方图。"""
    both = np.isfinite(py) & np.isfinite(mat) & (py != 0) & (mat != 0)
    diff = (py - mat)[both]
    if diff.size == 0:
        return

    fig, ax = plt.subplots(figsize=(10, 5))
    # 限制显示范围 ±99 分位
    p99 = np.percentile(np.abs(diff), 99)
    clip = max(0.5, p99 * 1.2)
    ax.hist(diff, bins=120, range=(-clip, clip),
            color='#4472C4', edgecolor='black', alpha=0.85)
    ax.axvline(0, color='k', lw=1.0)
    ax.axvline(diff.mean(), color='r', ls='--', lw=1.2,
               label=f'mean = {diff.mean():.4f} m')
    ax.axvline(diff.mean() + diff.std(), color='orange', ls=':', lw=1.0,
               label=f'mean ± std (std = {diff.std():.4f} m)')
    ax.axvline(diff.mean() - diff.std(), color='orange', ls=':', lw=1.0)
    ax.set_xlabel('Python − MATLAB  [m]')
    ax.set_ylabel('像元数')
    ax.set_title(f'差值分布直方图 — {component_label}', fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  已保存：{output_path}')


# -----------------------------------------------------------------------------
# 指标输出
# -----------------------------------------------------------------------------
def write_metrics_text(all_metrics: dict, output_path: str,
                       python_dir: str, matlab_dir: str,
                       py_meta: dict, mat_meta: dict):
    """汇总文本报告。"""
    lines = []
    lines.append('=' * 78)
    lines.append('Python vs MATLAB SM-VCE 结果对比报告')
    lines.append(f'生成时间：{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')
    lines.append('=' * 78)
    lines.append('')
    lines.append('数据源：')
    lines.append(f'  Python 目录：{python_dir}')
    lines.append(f'  MATLAB 目录：{matlab_dir}')
    lines.append('')
    lines.append('栅格元信息：')
    lines.append(f'  Python : shape = {py_meta["shape"]}, '
                 f'corner = ({py_meta["transform"].c:.6f}, {py_meta["transform"].f:.6f}), '
                 f'post = ({py_meta["transform"].a:.7f}, {py_meta["transform"].e:.7f})')
    lines.append(f'  MATLAB : shape = {mat_meta["shape"]}, '
                 f'corner = ({mat_meta["transform"].c:.6f}, {mat_meta["transform"].f:.6f}), '
                 f'post = ({mat_meta["transform"].a:.7f}, {mat_meta["transform"].e:.7f})')
    lines.append('')

    for key, label in COMPONENTS:
        m = all_metrics[key]
        lines.append('-' * 78)
        lines.append(f'分量：{label}  (enu_{key})')
        lines.append('-' * 78)
        lines.append(f'  覆盖统计：Python 有效 = {m["n_python_valid"]:>8d}，'
                     f'MATLAB 有效 = {m["n_matlab_valid"]:>8d}，'
                     f'共同有效 = {m["n_both_valid"]:>8d}  '
                     f'(占全图 {m["coverage_pct"]:.2f}%)')
        lines.append('')
        lines.append(f'  Python 均值 = {m["py_mean"]:+.5f} m，std = {m["py_std"]:.5f} m')
        lines.append(f'  MATLAB 均值 = {m["mat_mean"]:+.5f} m，std = {m["mat_std"]:.5f} m')
        lines.append('')
        lines.append(f'  差值统计 (Python − MATLAB)：')
        lines.append(f'    mean     = {m["mean_diff"]:+.5f} m')
        lines.append(f'    std      = {m["std_diff"]:.5f} m')
        lines.append(f'    RMSE     = {m["rmse"]:.5f} m')
        lines.append(f'    MAE      = {m["mae"]:.5f} m')
        lines.append(f'    max|d|   = {m["max_abs_diff"]:.5f} m')
        lines.append(f'    Pearson r= {m["corr"]:.5f}')
        lines.append(f'    |d|<5cm  = {m["pct_within_5cm"]:6.2f}%   '
                     f'|d|<10cm = {m["pct_within_10cm"]:6.2f}%   '
                     f'|d|<20cm = {m["pct_within_20cm"]:6.2f}%')
        lines.append('')

    lines.append('=' * 78)
    lines.append('指标说明：')
    lines.append('  mean_diff   差值均值 — 系统性偏差；接近 0 表示无整体偏移')
    lines.append('  std_diff    差值标准差 — 随机离散程度')
    lines.append('  RMSE        均方根误差 — 综合精度指标')
    lines.append('  MAE         平均绝对差 — 鲁棒精度指标')
    lines.append('  Pearson r   线性相关 — 1.0 = 完全线性一致')
    lines.append('  |d|<X       差值小于 X 米的像元占比 — 工程意义上的"基本相同"率')
    lines.append('=' * 78)

    with open(output_path, 'w', encoding='utf-8') as fp:
        fp.write('\n'.join(lines) + '\n')
    print(f'  已保存：{output_path}')


def write_metrics_csv(all_metrics: dict, output_path: str):
    fields = ['component', 'n_total', 'n_python_valid', 'n_matlab_valid', 'n_both_valid',
              'coverage_pct', 'py_mean', 'py_std', 'mat_mean', 'mat_std',
              'mean_diff', 'std_diff', 'rmse', 'mae', 'max_abs_diff', 'corr',
              'pct_within_5cm', 'pct_within_10cm', 'pct_within_20cm']
    with open(output_path, 'w', newline='', encoding='utf-8') as fp:
        w = csv.DictWriter(fp, fieldnames=fields)
        w.writeheader()
        for key, _ in COMPONENTS:
            row = {'component': key}
            row.update({k: all_metrics[key].get(k, '') for k in fields[1:]})
            w.writerow(row)
    print(f'  已保存：{output_path}')


def plot_summary(layers_py: dict, layers_mat: dict, transform,
                 all_metrics: dict, output_path: str,
                 vmin=-1.5, vmax=1.5, diff_lim=0.5):
    """总览图：3 行（E/N/U） × 3 列（Python / MATLAB / Diff）。"""
    row, col = layers_py['east'].shape
    extent = [transform.c, transform.c + transform.a * col,
              transform.f + transform.e * row, transform.f]

    fig, axes = plt.subplots(3, 3, figsize=(16, 18))
    norm_data = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
    norm_diff = TwoSlopeNorm(vmin=-diff_lim, vcenter=0, vmax=diff_lim)

    for ri, (key, label) in enumerate(COMPONENTS):
        py = layers_py[key]
        mat = layers_mat[key]
        diff = py - mat
        m = all_metrics[key]
        triples = [
            (py, f'Python {label}', norm_data, 'RdBu_r'),
            (mat, f'MATLAB {label}', norm_data, 'RdBu_r'),
            (diff, f'Diff  RMSE={m["rmse"]:.3f}m  r={m["corr"]:.3f}',
             norm_diff, 'RdBu_r'),
        ]
        for ci, (d, title, norm, cmap) in enumerate(triples):
            ax = axes[ri, ci]
            d_show = np.where(np.isfinite(d) & (d != 0), d, np.nan)
            im = ax.imshow(d_show, extent=extent, cmap=cmap, norm=norm,
                           aspect='auto', origin='upper', interpolation='nearest')
            ax.set_title(title, fontsize=11, fontweight='bold')
            if ci == 0:
                ax.set_ylabel(f'{label}\nLatitude (°N)')
            if ri == 2:
                ax.set_xlabel('Longitude (°E)')
            plt.colorbar(im, ax=ax, shrink=0.7, label='[m]')

    fig.suptitle('Python vs MATLAB SM-VCE — 三维形变场全景对比',
                 fontsize=15, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  已保存：{output_path}')


# -----------------------------------------------------------------------------
# main
# -----------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description='Python vs MATLAB SM-VCE 结果对比')
    parser.add_argument('--python-dir', type=str, default='',
                        help='Python 版 geotiff 目录（默认自动选最新）')
    parser.add_argument('--matlab-dir', type=str, default='',
                        help='MATLAB 版 geotiff 目录（默认自动选最新）')
    parser.add_argument('--out', type=str, default='',
                        help='输出目录（默认 Comparison_Python_vs_MATLAB_<ts>/）')
    parser.add_argument('--vmin', type=float, default=-1.5)
    parser.add_argument('--vmax', type=float, default=1.5)
    parser.add_argument('--diff-lim', type=float, default=0.5,
                        help='差值图的色标对称上限 (m)，默认 0.5')
    args = parser.parse_args()

    py_dir = args.python_dir or auto_locate_python_dir()
    mat_dir = args.matlab_dir or auto_locate_matlab_dir()
    print(f'Python 目录：{py_dir}')
    print(f'MATLAB 目录：{mat_dir}')

    if args.out:
        out_dir = args.out
    else:
        ts = datetime.now().strftime('%Y%m%d%H%M%S')
        out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               f'Comparison_Python_vs_MATLAB_{ts}')
    os.makedirs(out_dir, exist_ok=True)
    print(f'输出目录：{out_dir}\n')

    # 读取 3 个分量
    layers_py = {}
    layers_mat = {}
    py_meta = {}
    mat_meta = {}
    for key, _ in COMPONENTS:
        py_arr, py_tf, py_sh = load_layer(py_dir, key)
        mat_arr, mat_tf, mat_sh = load_layer(mat_dir, key)
        if py_sh != mat_sh:
            print(f'警告：{key} 分量栅格尺寸不一致 Python {py_sh} vs MATLAB {mat_sh}')
        layers_py[key] = py_arr
        layers_mat[key] = mat_arr
        py_meta = {'shape': py_sh, 'transform': py_tf}
        mat_meta = {'shape': mat_sh, 'transform': mat_tf}

    # 计算各分量指标
    print('正在计算指标 ...')
    all_metrics = {}
    for key, label in COMPONENTS:
        all_metrics[key] = compute_metrics(layers_py[key], layers_mat[key])
        m = all_metrics[key]
        print(f'  {label}: RMSE = {m["rmse"]:.4f} m，MAE = {m["mae"]:.4f} m，'
              f'r = {m["corr"]:.4f}，|d|<5cm = {m["pct_within_5cm"]:.1f}%')

    # 可视化
    print('\n正在绘图 ...')
    for key, label in COMPONENTS:
        plot_comparison(layers_py[key], layers_mat[key], py_meta['transform'],
                        key, label,
                        os.path.join(out_dir, f'Comparison_{key}.png'),
                        vmin=args.vmin, vmax=args.vmax, diff_lim=args.diff_lim)
        plot_scatter(layers_py[key], layers_mat[key], label,
                     os.path.join(out_dir, f'Scatter_{key}.png'))
        plot_diff_histogram(layers_py[key], layers_mat[key], label,
                            os.path.join(out_dir, f'Histogram_diff_{key}.png'))

    plot_summary(layers_py, layers_mat, py_meta['transform'], all_metrics,
                 os.path.join(out_dir, 'Summary.png'),
                 vmin=args.vmin, vmax=args.vmax, diff_lim=args.diff_lim)

    write_metrics_text(all_metrics, os.path.join(out_dir, 'metrics.txt'),
                       py_dir, mat_dir, py_meta, mat_meta)
    write_metrics_csv(all_metrics, os.path.join(out_dir, 'metrics.csv'))

    print(f'\n对比完成。结果目录：{out_dir}')


if __name__ == '__main__':
    main()
