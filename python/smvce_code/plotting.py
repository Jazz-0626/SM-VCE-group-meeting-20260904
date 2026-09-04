"""
Plotting utilities for SM-VCE results (publication-quality).
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from mpl_toolkits.axes_grid1 import make_axes_locatable
from .utils import lonlat2sub


# --------------- color map ------------------------------------------------
def _get_default_cmap():
    """Blue-White-Red diverging colormap suitable for displacement maps."""
    from matplotlib.colors import LinearSegmentedColormap
    colors = [
        (0.0, '#2166ac'),
        (0.25, '#67a9cf'),
        (0.5, '#f7f7f7'),
        (0.75, '#ef8a62'),
        (1.0, '#b2182b'),
    ]
    return LinearSegmentedColormap.from_list(
        'BlueWhiteRed', [(v, c) for v, c in colors], N=256)


def apply_publication_style():
    """Apply a restrained journal-style Matplotlib theme."""
    plt.rcParams.update({
        'font.family': 'DejaVu Serif',
        'mathtext.fontset': 'stix',
        'axes.unicode_minus': False,
        'axes.linewidth': 0.9,
        'xtick.direction': 'in',
        'ytick.direction': 'in',
        'xtick.major.width': 0.9,
        'ytick.major.width': 0.9,
        'xtick.major.size': 4.0,
        'ytick.major.size': 4.0,
        'savefig.dpi': 300,
    })


# --------------- main figure function -------------------------------------
def getfig(data, dls=None, flag_clb=True, flag_tklb=True,
           siz=None, lgdstr=None, cmap=None):
    """
    Create publication-quality multi-panel figure.

    Returns
    -------
    fig, axs, cbs
    """
    if data.ndim == 2:
        data = data[:, :, np.newaxis]
    row, col, figNum = data.shape

    # subplot layout
    if siz is None:
        nr = max(1, int(np.floor(np.sqrt(figNum))))
        nc = int(np.ceil(figNum / nr))
    else:
        nr, nc = siz

    # data limits
    if dls is None:
        dls = _auto_dls(data, figNum)
    else:
        dls = np.atleast_2d(dls)
        if dls.shape[0] == 1 and figNum > 1:
            dls = np.tile(dls, (figNum, 1))

    if cmap is None:
        cmap = _get_default_cmap()

    # aspect ratio from data
    aspect_ratio = row / col
    panel_w = 3.2  # inches per panel
    panel_h = panel_w * aspect_ratio
    fig_w = panel_w * nc + 1.0  # extra for colorbar / spacing
    fig_h = panel_h * nr + 0.8

    fig, axes = plt.subplots(nr, nc, figsize=(fig_w, fig_h), squeeze=False)
    axs = []
    cbs = []

    for idx in range(figNum):
        ri, ci_idx = divmod(idx, nc)
        ax = axes[ri, ci_idx]
        axs.append(ax)

        masked = np.ma.masked_invalid(data[:, :, idx])
        im = ax.imshow(masked, vmin=dls[idx, 0], vmax=dls[idx, 1],
                       cmap=cmap, aspect='equal', interpolation='nearest')

        # title
        if lgdstr is not None and idx < len(lgdstr):
            ax.set_title(lgdstr[idx], fontsize=10.5, fontweight='semibold', pad=6)

        # colorbar - use divider so it matches image height exactly
        if flag_clb:
            divider = make_axes_locatable(ax)
            cax = divider.append_axes('right', size='4%', pad=0.06)
            cb = fig.colorbar(im, cax=cax)
            cb.ax.tick_params(labelsize=8)
            cbs.append(cb)

        if not flag_tklb:
            ax.set_xticklabels([])
            ax.set_yticklabels([])

        ax.tick_params(axis='both', which='both', direction='in',
                       top=True, right=True, labelsize=8)

    # hide unused
    for idx in range(figNum, nr * nc):
        ri, ci_idx = divmod(idx, nc)
        axes[ri, ci_idx].set_visible(False)

    fig.subplots_adjust(wspace=0.35, hspace=0.25)
    return fig, axs, cbs


def compute_dls(data, symmetric=False):
    """Compute per-panel display limits from actual data values."""
    if data.ndim == 2:
        data = data[:, :, np.newaxis]
    return _compute_dls_from_list([data], symmetric=symmetric)


def compute_shared_dls(data_list, symmetric=False):
    """Compute per-panel display limits shared across multiple cubes."""
    cubes = []
    for data in data_list:
        if data.ndim == 2:
            data = data[:, :, np.newaxis]
        cubes.append(data)
    return _compute_dls_from_list(cubes, symmetric=symmetric)


# --------------- coordinate ticks -----------------------------------------
def coortick(ax, coor, dxy=(0.2, 0.2)):
    """Set lon/lat tick labels with degree symbols."""
    lon0 = coor['corner_lon']
    lat0 = coor['corner_lat']
    dlon = coor['post_lon']
    dlat = coor['post_lat']  # negative for north-up
    w = coor['width']
    h = coor['nlines']

    lon_min, lon_max = lon0, lon0 + dlon * (w - 1)
    lat_top, lat_bot = lat0, lat0 + dlat * (h - 1)
    lat_min, lat_max = min(lat_top, lat_bot), max(lat_top, lat_bot)

    # X ticks (longitude)
    Xs = np.ceil(lon_min / dxy[0]) * dxy[0]
    XX = np.arange(Xs, lon_max + dxy[0] * 0.5, dxy[0])
    XX = XX[(XX >= lon_min - 1e-9) & (XX <= lon_max + 1e-9)]
    xtick = (XX - lon0) / dlon

    # Y ticks (latitude)
    Ys = np.ceil(lat_min / dxy[1]) * dxy[1]
    YY = np.arange(Ys, lat_max + dxy[1] * 0.5, dxy[1])
    YY = YY[(YY >= lat_min - 1e-9) & (YY <= lat_max + 1e-9)]
    ytick = (YY - lat0) / dlat

    # labels
    xtl = [f'{x:.1f}°E' if x >= 0 else f'{-x:.1f}°W' for x in XX]
    ytl = [f'{y:.1f}°N' if y >= 0 else f'{-y:.1f}°S' for y in YY]

    # latitude labels should be ordered top→bottom matching pixel coords
    if dlat < 0:
        # top of image = lat0 (large lat), bottom = lat0+dlat*(h-1) (small lat)
        sort_idx = np.argsort(ytick)
        ytick = ytick[sort_idx]
        ytl = [ytl[i] for i in sort_idx]

    ax.set_xticks(xtick)
    ax.set_xticklabels(xtl, fontsize=8)
    ax.set_yticks(ytick)
    ax.set_yticklabels(ytl, fontsize=8, rotation=90, va='center')


# --------------- fault trace ----------------------------------------------
def plot_fault(ax, fault, coor, color='m', lw=1.5):
    """Plot fault traces."""
    if isinstance(fault, list):
        for fi in fault:
            subi = lonlat2sub(fi, coor)
            ax.plot(subi[:, 1], subi[:, 0], '-', color=color, linewidth=lw)


# --------------- 2D vector arrows -----------------------------------------
def select_sparse_vectors(enu, max_count=20, min_count=8,
                          candidate_quantile=15, relative_floor=0.08):
    """
    Select a sparse but spatially balanced set of horizontal displacement vectors.

    Parameters
    ----------
    enu : ndarray, shape (row, col, >=2)
        ENU displacement cube.
    max_count : int
        Maximum number of vectors kept for the whole panel.
    min_count : int
        Minimum number of vectors to keep when valid candidates exist.
    candidate_quantile : float
        Quantile used to suppress weak cell-wise candidates.
    relative_floor : float
        Additional relative threshold with respect to the strongest candidate.

    Returns
    -------
    ndarray, shape (N, 4)
        Sparse vectors in [col, row, east, north] format.
    """
    if enu.ndim != 3 or enu.shape[2] < 2 or max_count <= 0:
        return np.zeros((0, 4), dtype=float)

    east = np.asarray(enu[:, :, 0], dtype=float)
    north = np.asarray(enu[:, :, 1], dtype=float)
    mag = np.hypot(east, north)
    valid = np.isfinite(east) & np.isfinite(north) & np.isfinite(mag) & (mag > 0)
    if not np.any(valid):
        return np.zeros((0, 4), dtype=float)

    row, col = mag.shape
    grid_rows, grid_cols = _target_grid_shape(row, col, max_count)
    row_edges = np.linspace(0, row, grid_rows + 1, dtype=int)
    col_edges = np.linspace(0, col, grid_cols + 1, dtype=int)

    candidates = []
    for ir in range(grid_rows):
        for ic in range(grid_cols):
            r0, r1 = row_edges[ir], row_edges[ir + 1]
            c0, c1 = col_edges[ic], col_edges[ic + 1]
            if r1 <= r0 or c1 <= c0:
                continue

            valid_cell = valid[r0:r1, c0:c1]
            if not np.any(valid_cell):
                continue

            mag_cell = np.where(valid_cell, mag[r0:r1, c0:c1], -np.inf)
            rr, cc = np.unravel_index(np.argmax(mag_cell), mag_cell.shape)
            rr += r0
            cc += c0
            candidates.append([cc, rr, east[rr, cc], north[rr, cc], mag[rr, cc]])

    if not candidates:
        return np.zeros((0, 4), dtype=float)

    candidates = np.asarray(candidates, dtype=float)
    candidates = candidates[np.argsort(candidates[:, 4])[::-1]]

    mag_floor = max(
        np.percentile(candidates[:, 4], candidate_quantile),
        relative_floor * candidates[0, 4],
    )
    filtered = candidates[candidates[:, 4] >= mag_floor]

    keep_count = min(max_count, filtered.shape[0])
    if keep_count < min(min_count, candidates.shape[0]):
        keep_count = min(max_count, max(min_count, keep_count))
        filtered = candidates[:keep_count]
    else:
        filtered = filtered[:keep_count]

    return filtered[:, :4]


def plot2d_vector(ax, data, field_shape=None, color='#202020',
                  line_width=1,
                  mutation_scale=4,
                  alpha=0.9,
                  target_max_length_ratio=0.10,
                  add_reference=True):
    """
    Plot publication-style 2D displacement vectors.

    All arrows have the same head size and line width for a clean look.
    Arrow length encodes displacement magnitude; direction encodes azimuth.

    Parameters
    ----------
    data : ndarray (N, 4) - [col, row, east, north]
    field_shape : tuple(int, int), optional
        Raster shape as (row, col). Used to scale arrows and place the legend.
    """
    data = np.asarray(data, dtype=float)
    if data.size == 0:
        return

    mags = np.hypot(data[:, 2], data[:, 3])
    valid = np.isfinite(mags) & (mags > 0) & np.isfinite(data).all(axis=1)
    if not np.any(valid):
        return

    data = data[valid]
    mags = mags[valid]
    if data.shape[0] == 0:
        return

    if field_shape is None:
        row = max(1, int(np.nanmax(data[:, 1])) + 1)
        col = max(1, int(np.nanmax(data[:, 0])) + 1)
    else:
        row, col = field_shape
    panel_span = max(1.0, float(min(row, col)))

    mag_ref = np.percentile(mags, 95) if mags.size > 1 else mags[0]
    mag_ref = max(float(mag_ref), np.finfo(float).eps)
    scale = target_max_length_ratio * panel_span / mag_ref

    for i in range(data.shape[0]):
        xy0 = data[i, :2]
        dxy = data[i, 2:4]
        dxy_s = dxy * scale

        ann = ax.annotate(
            '',
            xy=(xy0[0] + dxy_s[0], xy0[1] - dxy_s[1]),
            xytext=(xy0[0], xy0[1]),
            arrowprops=dict(
                arrowstyle='-|>',
                color=color,
                lw=line_width,
                mutation_scale=mutation_scale,
                alpha=alpha,
                shrinkA=0,
                shrinkB=0,
                joinstyle='round',
                capstyle='round',
            ),
            zorder=5,
        )
        if ann.arrow_patch is not None:
            ann.arrow_patch.set_path_effects([
                pe.Stroke(linewidth=line_width + 1.0, foreground='white', alpha=0.6),
                pe.Normal(),
            ])

    if add_reference:
        ref_mag = _nice_number(np.percentile(mags, 75))
        if ref_mag > 0:
            _draw_reference_arrow(
                ax=ax,
                ref_mag=ref_mag,
                scale=scale,
                row=row,
                col=col,
                color=color,
                line_width=line_width,
                mutation_scale=mutation_scale,
            )


# --------------- colorbar label -------------------------------------------
def clbtitle(cb, title_str, fontsize=9):
    """Set colorbar label."""
    cb.set_label(title_str, fontsize=fontsize, rotation=90, labelpad=4)


# --------------- save -----------------------------------------------------
def save_figure(fig, filepath, dpi=300):
    """Save figure with tight bounding box."""
    fig.savefig(filepath, dpi=dpi, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    print(f'已保存：{filepath}')


# --------------- helpers --------------------------------------------------
def _target_grid_shape(row, col, max_count):
    """Choose a grid close to the image aspect ratio while keeping <= max_count cells."""
    aspect = row / max(col, 1)
    grid_rows = max(1, int(round(np.sqrt(max_count * aspect))))
    grid_cols = max(1, int(np.ceil(max_count / grid_rows)))
    while grid_rows * grid_cols > max_count:
        if grid_cols >= grid_rows and grid_cols > 1:
            grid_cols -= 1
        elif grid_rows > 1:
            grid_rows -= 1
        else:
            break
    while grid_rows * grid_cols < max_count:
        if grid_cols <= grid_rows:
            if grid_rows * (grid_cols + 1) <= max_count:
                grid_cols += 1
            else:
                break
        else:
            if (grid_rows + 1) * grid_cols <= max_count:
                grid_rows += 1
            else:
                break
    return grid_rows, grid_cols


def _compute_dls_from_list(cubes, symmetric=False):
    fig_num = cubes[0].shape[2]
    dls = []
    for idx in range(fig_num):
        valid_parts = []
        for cube in cubes:
            d = cube[:, :, idx]
            valid = d[np.isfinite(d)]
            if valid.size > 0:
                valid_parts.append(valid)

        if not valid_parts:
            dls.append([-0.5, 0.5] if symmetric else [0.0, 1.0])
            continue

        valid = np.concatenate(valid_parts)
        if symmetric:
            vmax = float(np.nanmax(np.abs(valid)))
            if vmax == 0:
                vmax = 0.5
            dls.append([-vmax, vmax])
        else:
            vmin = float(np.nanmin(valid))
            vmax = float(np.nanmax(valid))
            if vmin == vmax:
                pad = 0.5 if vmin == 0 else abs(vmin) * 0.05
                vmin -= pad
                vmax += pad
            dls.append([vmin, vmax])
    return np.array(dls)


def _auto_dls(data, figNum):
    dls = []
    for i in range(figNum):
        d = data[:, :, i]
        vmin = np.nanmin(d) if np.any(~np.isnan(d)) else -0.5
        vmax = np.nanmax(d) if np.any(~np.isnan(d)) else 0.5
        if vmin == vmax:
            vmin -= 0.5
            vmax += 0.5
        dls.append([vmin, vmax])
    return np.array(dls)


def _nice_number(value):
    """Round a positive value to a compact legend-friendly number."""
    if not np.isfinite(value) or value <= 0:
        return 0.0
    exponent = np.floor(np.log10(value))
    fraction = value / (10 ** exponent)
    if fraction < 1.5:
        nice_fraction = 1.0
    elif fraction < 3.5:
        nice_fraction = 2.0
    elif fraction < 7.5:
        nice_fraction = 5.0
    else:
        nice_fraction = 10.0
    return nice_fraction * (10 ** exponent)


def _draw_reference_arrow(ax, ref_mag, scale, row, col, color, line_width, mutation_scale):
    """Draw a small horizontal reference arrow in the lower-left corner."""
    x0 = 0.07 * (col - 1)
    y0 = 0.90 * (row - 1)
    dx = ref_mag * scale

    ann = ax.annotate(
        '',
        xy=(x0 + dx, y0),
        xytext=(x0, y0),
        arrowprops=dict(
            arrowstyle='-|>',
            color=color,
            lw=line_width,
            mutation_scale=mutation_scale,
            alpha=0.95,
            shrinkA=0,
            shrinkB=0,
            joinstyle='round',
            capstyle='round',
        ),
        zorder=6,
    )
    if ann.arrow_patch is not None:
        ann.arrow_patch.set_path_effects([
            pe.Stroke(linewidth=line_width + 1.0, foreground='white', alpha=0.75),
            pe.Normal(),
        ])

    ax.text(
        x0,
        y0 + 0.045 * row,
        f'{ref_mag:.2g} m',
        fontsize=8,
        color=color,
        ha='left',
        va='bottom',
        zorder=6,
        bbox=dict(boxstyle='round,pad=0.18', facecolor='white', edgecolor='none', alpha=0.72),
    )
