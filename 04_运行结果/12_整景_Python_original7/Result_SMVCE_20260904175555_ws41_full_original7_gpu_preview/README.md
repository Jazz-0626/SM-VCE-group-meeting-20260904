# SM-VCE 结果说明

- 生成时间：`2026-09-04 17:56:04`
- 结果根目录：`E:\#Post\2026.09.04 今日组会\04_运行结果\12_整景_Python_original7\Result_SMVCE_20260904175555_ws41_full_original7_gpu_preview`
- GeoTIFF 目录：`E:\#Post\2026.09.04 今日组会\04_运行结果\12_整景_Python_original7\Result_SMVCE_20260904175555_ws41_full_original7_gpu_preview\geotiff`
- GRD 导出：`disabled`
- SM-VCE 求解耗时：`679.56 s`
- 栅格尺寸：`1027 x 1373`
- 坐标参考：`+proj=longlat +datum=WGS84 +no_defs`
- 左上角坐标：`(100.95363205000001, 37.965229900000004)`
- 像元间距：`dlon=0.0004329, dlat=-0.0003882`

## 数据含义

- `enu_*`：SM-VCE 最终解算得到的东、北、天三维形变，单位为 `m`。
- `enu_std_*`：对应 ENU 三个分量的中误差，单位为 `m`。
- `obs_std_*`：每组原始观测在局部 VCE 估计后的观测中误差，单位为 `m`。
- `sita_*`：每组观测对应的方差分量估计值。
- `shpcount_*`：每个局部窗口内参与求解的有效同质样本数。
- `para_sm_*`：与位移一起估计得到的局部应变模型参数。

## 输入观测对应关系

| 序号 | 输入文件 | 几何标记 |
| --- | --- | --- |
| 1 | `S1_As_t128_DInSAR.disp.tif` | `1` |
| 2 | `S1_As_t128_POT_LOS.disp.tif` | `1` |
| 3 | `S1_Des_DInSAR.disp.tif` | `1` |
| 4 | `S1_Des_POT_LOS.disp.tif` | `1` |
| 5 | `ALOS2_DInSAR_LOS.disp.tif` | `1` |
| 6 | `ALOS2_Des_POT_AZI.disp.tif` | `2` |
| 7 | `ALOS2_Des_MAI.disp.tif` | `2` |

## 输出图层清单

| 变量名 | GeoTIFF | GRD | 单位 | 说明 |
| --- | --- | --- | --- | --- |
| `enu_east` | `enu_east.tif` | `未生成` | `m` | 东向形变 |
| `enu_north` | `enu_north.tif` | `未生成` | `m` | 北向形变 |
| `enu_up` | `enu_up.tif` | `未生成` | `m` | 垂向形变 |
| `enu_std_east` | `enu_std_east.tif` | `未生成` | `m` | 东向形变中误差 |
| `enu_std_north` | `enu_std_north.tif` | `未生成` | `m` | 北向形变中误差 |
| `enu_std_up` | `enu_std_up.tif` | `未生成` | `m` | 垂向形变中误差 |
| `para_sm_e_dx` | `para_sm_e_dx.tif` | `未生成` | `model_parameter` | 应变模型参数：para_sm_e_dx |
| `para_sm_e_dy` | `para_sm_e_dy.tif` | `未生成` | `model_parameter` | 应变模型参数：para_sm_e_dy |
| `para_sm_n_dx` | `para_sm_n_dx.tif` | `未生成` | `model_parameter` | 应变模型参数：para_sm_n_dx |
| `para_sm_n_dy` | `para_sm_n_dy.tif` | `未生成` | `model_parameter` | 应变模型参数：para_sm_n_dy |
| `para_sm_u_dx` | `para_sm_u_dx.tif` | `未生成` | `model_parameter` | 应变模型参数：para_sm_u_dx |
| `para_sm_u_dy` | `para_sm_u_dy.tif` | `未生成` | `model_parameter` | 应变模型参数：para_sm_u_dy |
| `obs_std_01_S1_As_t128_DInSAR.disp` | `obs_std_01_S1_As_t128_DInSAR.disp.tif` | `未生成` | `m` | 第 1 组观测中误差：S1_As_t128_DInSAR.disp.tif |
| `sita_01_S1_As_t128_DInSAR.disp` | `sita_01_S1_As_t128_DInSAR.disp.tif` | `未生成` | `unit_weight_variance` | 第 1 组观测的方差分量：S1_As_t128_DInSAR.disp.tif |
| `shpcount_01_S1_As_t128_DInSAR.disp` | `shpcount_01_S1_As_t128_DInSAR.disp.tif` | `未生成` | `count` | 第 1 组观测的有效同质样本数：S1_As_t128_DInSAR.disp.tif |
| `obs_std_02_S1_As_t128_POT_LOS.disp` | `obs_std_02_S1_As_t128_POT_LOS.disp.tif` | `未生成` | `m` | 第 2 组观测中误差：S1_As_t128_POT_LOS.disp.tif |
| `sita_02_S1_As_t128_POT_LOS.disp` | `sita_02_S1_As_t128_POT_LOS.disp.tif` | `未生成` | `unit_weight_variance` | 第 2 组观测的方差分量：S1_As_t128_POT_LOS.disp.tif |
| `shpcount_02_S1_As_t128_POT_LOS.disp` | `shpcount_02_S1_As_t128_POT_LOS.disp.tif` | `未生成` | `count` | 第 2 组观测的有效同质样本数：S1_As_t128_POT_LOS.disp.tif |
| `obs_std_03_S1_Des_DInSAR.disp` | `obs_std_03_S1_Des_DInSAR.disp.tif` | `未生成` | `m` | 第 3 组观测中误差：S1_Des_DInSAR.disp.tif |
| `sita_03_S1_Des_DInSAR.disp` | `sita_03_S1_Des_DInSAR.disp.tif` | `未生成` | `unit_weight_variance` | 第 3 组观测的方差分量：S1_Des_DInSAR.disp.tif |
| `shpcount_03_S1_Des_DInSAR.disp` | `shpcount_03_S1_Des_DInSAR.disp.tif` | `未生成` | `count` | 第 3 组观测的有效同质样本数：S1_Des_DInSAR.disp.tif |
| `obs_std_04_S1_Des_POT_LOS.disp` | `obs_std_04_S1_Des_POT_LOS.disp.tif` | `未生成` | `m` | 第 4 组观测中误差：S1_Des_POT_LOS.disp.tif |
| `sita_04_S1_Des_POT_LOS.disp` | `sita_04_S1_Des_POT_LOS.disp.tif` | `未生成` | `unit_weight_variance` | 第 4 组观测的方差分量：S1_Des_POT_LOS.disp.tif |
| `shpcount_04_S1_Des_POT_LOS.disp` | `shpcount_04_S1_Des_POT_LOS.disp.tif` | `未生成` | `count` | 第 4 组观测的有效同质样本数：S1_Des_POT_LOS.disp.tif |
| `obs_std_05_ALOS2_DInSAR_LOS.disp` | `obs_std_05_ALOS2_DInSAR_LOS.disp.tif` | `未生成` | `m` | 第 5 组观测中误差：ALOS2_DInSAR_LOS.disp.tif |
| `sita_05_ALOS2_DInSAR_LOS.disp` | `sita_05_ALOS2_DInSAR_LOS.disp.tif` | `未生成` | `unit_weight_variance` | 第 5 组观测的方差分量：ALOS2_DInSAR_LOS.disp.tif |
| `shpcount_05_ALOS2_DInSAR_LOS.disp` | `shpcount_05_ALOS2_DInSAR_LOS.disp.tif` | `未生成` | `count` | 第 5 组观测的有效同质样本数：ALOS2_DInSAR_LOS.disp.tif |
| `obs_std_06_ALOS2_Des_POT_AZI.disp` | `obs_std_06_ALOS2_Des_POT_AZI.disp.tif` | `未生成` | `m` | 第 6 组观测中误差：ALOS2_Des_POT_AZI.disp.tif |
| `sita_06_ALOS2_Des_POT_AZI.disp` | `sita_06_ALOS2_Des_POT_AZI.disp.tif` | `未生成` | `unit_weight_variance` | 第 6 组观测的方差分量：ALOS2_Des_POT_AZI.disp.tif |
| `shpcount_06_ALOS2_Des_POT_AZI.disp` | `shpcount_06_ALOS2_Des_POT_AZI.disp.tif` | `未生成` | `count` | 第 6 组观测的有效同质样本数：ALOS2_Des_POT_AZI.disp.tif |
| `obs_std_07_ALOS2_Des_MAI.disp` | `obs_std_07_ALOS2_Des_MAI.disp.tif` | `未生成` | `m` | 第 7 组观测中误差：ALOS2_Des_MAI.disp.tif |
| `sita_07_ALOS2_Des_MAI.disp` | `sita_07_ALOS2_Des_MAI.disp.tif` | `未生成` | `unit_weight_variance` | 第 7 组观测的方差分量：ALOS2_Des_MAI.disp.tif |
| `shpcount_07_ALOS2_Des_MAI.disp` | `shpcount_07_ALOS2_Des_MAI.disp.tif` | `未生成` | `count` | 第 7 组观测的有效同质样本数：ALOS2_Des_MAI.disp.tif |

## 备注

- GeoTIFF 始终生成；GRD 仅在启用且系统安装 GMT 时生成。
- 程序内部虽然使用 `var` 这个字段名，但导出的 `enu_std_*` 与 `obs_std_*` 实际表示中误差，不是方差。
- 数据中的 `NaN` 表示该像元无效或该位置没有成功完成求解。
- 每个 GeoTIFF 或 GRD 文件只保存一个二维图层，便于单独调用和后处理。
