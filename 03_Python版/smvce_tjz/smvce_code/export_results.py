"""
Export SM-VCE results to GeoTIFF and GRD files with a human-readable manifest.
"""
import os
import re
from datetime import datetime

import numpy as np
import rasterio
import xarray as xr
from rasterio.transform import Affine


def export_smvce_results(result_smvce, coor, datainfo, output_root, write_grd=True):
    """
    Export SM-VCE outputs into GeoTIFF and GRD folders.

    Parameters
    ----------
    result_smvce : dict
        Output dictionary returned by smvce_solve3d.
    coor : dict
        Coordinate information.
    datainfo : list
        Parsed rows from SMVCE_DATA/data_information.
    output_root : str
        Root output folder. Two subfolders will be created: geotiff and grd.

    Returns
    -------
    dict
        Paths of created outputs.
    """
    geotiff_dir = os.path.join(output_root, 'geotiff')
    grd_dir = os.path.join(output_root, 'grd')
    os.makedirs(geotiff_dir, exist_ok=True)
    if write_grd:
        os.makedirs(grd_dir, exist_ok=True)

    exports = []
    enu = np.asarray(result_smvce['enu'], dtype=np.float64)
    enu_std = np.asarray(result_smvce['var']['enu'], dtype=np.float64)
    obs_std = np.asarray(result_smvce['var']['obs'], dtype=np.float64)
    para_sm = result_smvce['para_sm']
    sita = result_smvce['sita']
    shpcount = result_smvce['SHPcount']
    input_data = result_smvce.get('InputData', {})

    exports.extend([
        ('enu_east', enu[:, :, 0], 'm', '东向形变'),
        ('enu_north', enu[:, :, 1], 'm', '北向形变'),
        ('enu_up', enu[:, :, 2], 'm', '垂向形变'),
        ('enu_std_east', enu_std[:, :, 0], 'm', '东向形变中误差'),
        ('enu_std_north', enu_std[:, :, 1], 'm', '北向形变中误差'),
        ('enu_std_up', enu_std[:, :, 2], 'm', '垂向形变中误差'),
    ])

    para_names = _para_sm_names(
        para_sm.shape[2],
        input_data.get('fsmpara'),
        input_data.get('flag_if_2D')
    )
    for idx, para_name in enumerate(para_names):
        exports.append((
            para_name,
            para_sm[:, :, idx],
            'model_parameter',
            f'应变模型参数：{para_name}'
        ))

    for idx, info in enumerate(datainfo):
        obs_label = _safe_label(info[0])
        exports.extend([
            (
                f'obs_std_{idx + 1:02d}_{obs_label}',
                obs_std[:, :, idx],
                'm',
                f'第 {idx + 1} 组观测中误差：{info[0]}'
            ),
            (
                f'sita_{idx + 1:02d}_{obs_label}',
                sita[:, :, idx],
                'unit_weight_variance',
                f'第 {idx + 1} 组观测的方差分量：{info[0]}'
            ),
            (
                f'shpcount_{idx + 1:02d}_{obs_label}',
                shpcount[:, :, idx],
                'count',
                f'第 {idx + 1} 组观测的有效同质样本数：{info[0]}'
            ),
        ])

    written = []
    for var_name, data_2d, unit, description in exports:
        tif_path = os.path.join(geotiff_dir, f'{var_name}.tif')
        grd_path = os.path.join(grd_dir, f'{var_name}.grd')
        _write_geotiff(tif_path, data_2d, coor)
        grd_name = ''
        if write_grd:
            _write_grd(grd_path, data_2d, coor, var_name, unit, description)
            grd_name = os.path.basename(grd_path)
        written.append((var_name, os.path.basename(tif_path), grd_name, unit, description))

    readme_path = os.path.join(output_root, 'README.md')
    _write_readme(readme_path, output_root, geotiff_dir, grd_dir, datainfo, written, result_smvce, write_grd)

    return {
        'output_root': output_root,
        'geotiff_dir': geotiff_dir,
        'grd_dir': grd_dir,
        'readme_path': readme_path,
    }


def _write_geotiff(path, data_2d, coor):
    """Write a single-band GeoTIFF."""
    data_2d = np.asarray(data_2d, dtype=np.float32)
    height, width = data_2d.shape
    transform = Affine(
        coor['post_lon'], 0.0, coor['corner_lon'],
        0.0, coor['post_lat'], coor['corner_lat']
    )
    crs = '+proj=longlat +datum=WGS84 +no_defs'

    with rasterio.open(
        path,
        'w',
        driver='GTiff',
        height=height,
        width=width,
        count=1,
        dtype='float32',
        crs=crs,
        transform=transform,
        nodata=np.nan,
        compress='deflate',
        predictor=3,
        zlevel=4,
    ) as dst:
        dst.write(data_2d, 1)


def _write_grd(path, data_2d, coor, var_name, unit, description):
    """Write a GMT-compatible binary grid (.grd) via gmt xyz2grd."""
    import subprocess, tempfile
    data_2d = np.asarray(data_2d, dtype=np.float32)
    height, width = data_2d.shape

    lon_min = coor['corner_lon']
    lon_max = lon_min + (width - 1) * coor['post_lon']
    lat_max = coor['corner_lat']  # corner_lat is north edge
    lat_min = lat_max + (height - 1) * coor['post_lat']  # post_lat is negative
    inc_lon = abs(coor['post_lon'])
    inc_lat = abs(coor['post_lat'])

    # Flip data so rows go south→north (GMT convention)
    data_sn = data_2d[::-1, :].copy()
    # Replace NaN with GMT NaN convention
    data_sn = np.where(np.isfinite(data_sn), data_sn, np.nan)

    # Write binary float32 (row-major, south→north, west→east)
    tmp_bin = tempfile.mktemp(suffix='.bin')
    data_sn.astype(np.float32).tofile(tmp_bin)

    R = f'-R{lon_min}/{lon_max}/{lat_min}/{lat_max}'
    I = f'-I{inc_lon}/{inc_lat}'

    try:
        subprocess.run(
            f'gmt xyz2grd {tmp_bin} {R} {I} -G{path} -ZBLf -fg '
            f'-D+x"longitude"+y"latitude"+z"{description} [{unit}]"+t"{var_name}"',
            shell=True, check=True, capture_output=True
        )
    finally:
        if os.path.exists(tmp_bin):
            os.unlink(tmp_bin)


def _write_readme(path, output_root, geotiff_dir, grd_dir, datainfo, written, result_smvce, write_grd):
    """Write a concise manifest for all exported layers."""
    total_time = float(result_smvce.get('total_time', np.nan))
    now_str = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    coor = result_smvce.get('coor', {})
    enu = np.asarray(result_smvce.get('enu'))
    row = int(enu.shape[0]) if enu.ndim >= 2 else 0
    col = int(enu.shape[1]) if enu.ndim >= 2 else 0
    crs = '+proj=longlat +datum=WGS84 +no_defs'
    lines = [
        '# SM-VCE 结果说明',
        '',
        f'- 生成时间：`{now_str}`',
        f'- 结果根目录：`{output_root}`',
        f'- GeoTIFF 目录：`{geotiff_dir}`',
        f'- GRD 导出：`{"enabled" if write_grd else "disabled"}`',
        f'- SM-VCE 求解耗时：`{total_time:.2f} s`',
        f'- 栅格尺寸：`{row} x {col}`',
        f'- 坐标参考：`{crs}`',
        f'- 左上角坐标：`({coor.get("corner_lon")}, {coor.get("corner_lat")})`',
        f'- 像元间距：`dlon={coor.get("post_lon")}, dlat={coor.get("post_lat")}`',
        '',
        '## 数据含义',
        '',
        '- `enu_*`：SM-VCE 最终解算得到的东、北、天三维形变，单位为 `m`。',
        '- `enu_std_*`：对应 ENU 三个分量的中误差，单位为 `m`。',
        '- `obs_std_*`：每组原始观测在局部 VCE 估计后的观测中误差，单位为 `m`。',
        '- `sita_*`：每组观测对应的方差分量估计值。',
        '- `shpcount_*`：每个局部窗口内参与求解的有效同质样本数。',
        '- `para_sm_*`：与位移一起估计得到的局部应变模型参数。',
        '',
        '## 输入观测对应关系',
        '',
        '| 序号 | 输入文件 | 几何标记 |',
        '| --- | --- | --- |',
    ]
    for idx, info in enumerate(datainfo, start=1):
        geom = info[1] if len(info) > 1 else ''
        lines.append(f'| {idx} | `{info[0]}` | `{geom}` |')

    lines.extend([
        '',
        '## 输出图层清单',
        '',
        '| 变量名 | GeoTIFF | GRD | 单位 | 说明 |',
        '| --- | --- | --- | --- | --- |',
    ])
    for var_name, tif_name, grd_name, unit, description in written:
        lines.append(f'| `{var_name}` | `{tif_name}` | `{grd_name or "未生成"}` | `{unit}` | {description} |')

    lines.extend([
        '',
        '## 备注',
        '',
        '- GeoTIFF 始终生成；GRD 仅在启用且系统安装 GMT 时生成。',
        '- 程序内部虽然使用 `var` 这个字段名，但导出的 `enu_std_*` 与 `obs_std_*` 实际表示中误差，不是方差。',
        '- 数据中的 `NaN` 表示该像元无效或该位置没有成功完成求解。',
        '- 每个 GeoTIFF 或 GRD 文件只保存一个二维图层，便于单独调用和后处理。',
    ])

    with open(path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines) + '\n')


def _para_sm_names(count, fsmpara, flag_if_2D):
    """Return readable names for the exported strain-model parameters."""
    if fsmpara == 2 and flag_if_2D == 0 and count == 6:
        return [
            'para_sm_e_dx',
            'para_sm_e_dy',
            'para_sm_n_dx',
            'para_sm_n_dy',
            'para_sm_u_dx',
            'para_sm_u_dy',
        ]
    if fsmpara == 2 and flag_if_2D == 1 and count == 4:
        return [
            'para_sm_e_dx',
            'para_sm_e_dy',
            'para_sm_u_dx',
            'para_sm_u_dy',
        ]
    if fsmpara == 3 and flag_if_2D == 0 and count == 9:
        return [
            'para_sm_e_dx',
            'para_sm_e_dy',
            'para_sm_e_dz',
            'para_sm_n_dx',
            'para_sm_n_dy',
            'para_sm_n_dz',
            'para_sm_u_dx',
            'para_sm_u_dy',
            'para_sm_u_dz',
        ]
    if fsmpara == 3 and flag_if_2D == 1 and count == 6:
        return [
            'para_sm_e_dx',
            'para_sm_e_dy',
            'para_sm_e_dz',
            'para_sm_u_dx',
            'para_sm_u_dy',
            'para_sm_u_dz',
        ]
    return [f'para_sm_{idx + 1:02d}' for idx in range(count)]


def _safe_label(name):
    """Make a filename-safe label."""
    stem = os.path.splitext(os.path.basename(name))[0]
    return re.sub(r'[^A-Za-z0-9._-]+', '_', stem)
