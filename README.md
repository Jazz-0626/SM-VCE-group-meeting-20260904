# SM-VCE 组会资料与代码

本仓库整理了 2026 年 9 月 4 日组会中使用的 SM-VCE（Strain Model–Variance Component Estimation，应变模型—方差分量估计）文献、MATLAB 原版程序、Python 移植版、运行结果和最终汇报 PPT。

SM-VCE 用局部应变模型描述窗口内连续形变，并通过方差分量估计为不同来源、不同精度的 LOS/AZI 观测自适应定权，最终反演东西、南北和垂直三个方向的地表位移。

## 仓库内容

| 目录 | 内容 |
|---|---|
| `01_文献/` | 核心论文、参考文献表及文献获取说明 |
| `02_MATLAB原版/` | 刘计洪提供的 MATLAB 原版程序、工具手册和组会实验入口 |
| `03_Python版/` | Python 移植版、震例适配代码及结果对比工具 |
| `04_运行结果/` | MATLAB/Python 测试、整景结果、GeoTIFF、GRD、结果图和定量对比表 |
| `05_PPT/` | 最终版 33 页组会 PPT |

`06_讲稿与说明/` 为本地讲稿目录，不纳入公开仓库。

## 克隆与下载

论文 PDF、MAT、GeoTIFF、GRD 和 PNG 等二进制文件由 Git LFS 管理。首次使用前请安装 [Git LFS](https://git-lfs.com/)，否则克隆后可能只能看到指针文件。

```powershell
git lfs install
git clone https://github.com/Jazz-0626/SM-VCE-group-meeting-20260904.git
cd SM-VCE-group-meeting-20260904
git lfs pull
```

## 数据说明

受数据授权和体积限制，仓库不提供原始 InSAR/POT 观测，也不提供由原始观测重新打包得到的 `DATA*.mat`。仓库中发布的是文献和已经生成的派生结果。

因此，克隆仓库后可以直接查看代码、PPT 与既有结果，但若要重新运行震例，需要另行取得原始数据，并将其放入本地的 `SMVCE_DATA/` 目录。该目录和 `DATA*.mat` 已由 `.gitignore` 排除，请勿提交。

## Python 版快速使用

### 1. 创建环境

```powershell
cd '03_Python版\smvce_tjz'
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

程序可在 CPU 上运行；若希望使用 CUDA 快速路径，还需根据本机 CUDA 版本另行安装 PyTorch。

### 2. 查看参数并运行

```powershell
python smvce_main.py --help
python smvce_main.py `
  --data-dir 'D:\path\to\SMVCE_DATA' `
  --windowsize 41 `
  --output-dir '.\output'
```

无 CUDA 环境时增加 `--no-cuda`。建议先用 `--roi R1 R2 C1 C2` 做小范围测试，确认数据格式和参数后再运行整景。

> `--preliminary-gpu-no-fault-separation` 只用于快速生成整景预览。该模式不执行跨断层邻域剔除，不能视为严格参考解。

更完整的输入格式、算法公式和参数说明见 [`03_Python版/smvce_tjz/README.md`](03_Python版/smvce_tjz/README.md)。

## MATLAB 原版快速使用

1. 在 MATLAB 中进入 `02_MATLAB原版/SMVCE_MATLAB_原版/`。
2. 将 `SMVCE_code/` 添加到 MATLAB 路径。
3. 按工具手册准备本地 `SMVCE_DATA/` 数据目录。
4. 运行 `SMVCE_main.m`。

若要复现本次组会的批量实验，可使用 `02_MATLAB原版/run_groupmeeting_experiment.m`。该脚本将结果写入原版源码目录之外，避免修改刘计洪提供的程序。

## 两个版本的定位

- MATLAB 版是刘计洪提供的原始实现，作为算法基准和方法来源。
- Python 版是面向实际地震震例的移植与适配版本，加入了命令行接口、结果导出、可选 CUDA 加速及震例相关处理。
- Python 版的震例修正不是对所有数据都成立的通用默认参数。窗口大小、断层几何、缺失值分布、观测组合和条带噪声应针对每个震例重新评估。

## 运行结果

`04_运行结果/` 同时保留了：

- 7 组原始观测与加入高噪声观测后的 9 组对比；
- MATLAB 原版、Python 复现模式和 Python 震例特殊处理结果；
- ROI 测试和全尺寸整景结果；
- ENU 位移、标准差、方差分量、有效样本数、GeoTIFF/GRD 及定量对比图表。

其中带有 `gpu_preview` 的结果属于 GPU 快速预览，应结合相应 `run_config.json` 和结果说明解读。

## 文献与引用

核心文献的 DOI、开放获取状态和推荐阅读顺序见 [`01_文献/文献说明.md`](01_文献/文献说明.md)，BibTeX 条目见 [`01_文献/refs.bib`](01_文献/refs.bib)。仓库内两篇论文 PDF 均为开放获取版本；引用算法或论文时，请使用原始论文信息，而不是只引用本仓库。

## 使用与责任说明

本仓库用于科研交流和方法复现。MATLAB SM-VCE 程序归属及引用请遵循原作者和随附材料的说明；论文、数据和第三方组件仍适用各自的许可条件。使用结果开展科研分析时，请自行核验参数、收敛情况、空间覆盖范围及物理合理性。
