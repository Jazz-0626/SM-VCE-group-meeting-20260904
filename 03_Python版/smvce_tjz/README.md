# SM-VCE Python 工程详细说明

本工程实现了基于应变模型方差分量估计（Strain Model Variance Component Estimation, SM-VCE）的多源InSAR/POT观测三维地表形变反演方法。该方法从MATLAB版本转换而来，主要用于大地测量学和地震学领域的地表形变监测。

## 目录

1. [项目概述](#1-项目概述)
2. [算法原理与数学公式](#2-算法原理与数学公式)
3. [详细代码流程](#3-详细代码流程)
4. [核心模块详解](#4-核心模块详解)
5. [输入输出约定](#5-输入输出约定)
6. [使用方法](#6-使用方法)
7. [工程限制与注意事项](#7-工程限制与注意事项)

---

## 1. 项目概述

### 1.1 项目背景

在大地测量学中，InSAR（合成孔径雷达干涉测量）和POT（像素偏移跟踪）技术可以获取地表形变的观测数据。这些观测通常包含：
- **LOS观测**：视线方向（Line-of-Sight）形变分量
- **AZI观测**：方位向形变分量
- **E-W观测**：东西向形变分量
- **N-S观测**：南北向形变分量

由于不同观测类型的几何特性和误差特性不同，需要一种能够联合处理多源观测并自适应估计观测权重的方法。SM-VCE方法正是为解决这一问题而设计。

### 1.2 核心思想

SM-VCE方法的核心思想包括：

1. **应变模型（Strain Model, SM）**：假设局部区域内形变可以表示为位移加上应变参数的线性函数
2. **方差分量估计（Variance Component Estimation, VCE）**：通过迭代估计各类观测的方差分量，自适应调整观测权重
3. **窗口化处理**：在每个像元周围建立局部窗口，利用窗口内的观测进行约束求解

### 1.3 主要功能

- 支持多种观测类型（LOS、AZI、E-W、N-S、Vertical）的联合处理
- 经典WLS（加权最小二乘）三维反演
- SM-VCE自适应加权三维反演
- 支持二维（E-W + Vertical）和三维（E-W + N-S + Vertical）求解
- 支持二维和三维应变模型
- 同质点筛选（SMAD）
- 自适应窗口大小
- 断层约束处理
- CUDA加速（可选）

---

## 2. 算法原理与数学公式

### 2.1 观测方程建立

#### 2.1.1 观测几何系数

对于每种观测类型，建立其与三维位移（E, N, U）的线性关系：

**LOS观测**（视线方向）：

$$d_{LOS} = -\epsilon \cdot \sin(\theta) \cdot \sin(\alpha) \cdot E - \epsilon \cdot \sin(\theta) \cdot \cos(\alpha) \cdot N + \cos(\theta) \cdot U$$

其中：
- $\theta$：入射角（incidence angle）
- $\alpha$：方位角（azimuth angle），相对于北方向
- $\epsilon$：左右视标记（右视=1，左视=-1）

**AZI观测**（方位向）：

$$d_{AZI} = -\cos(\alpha) \cdot E + \sin(\alpha) \cdot N$$

**E-W观测**：

$$d_{EW} = E$$

**N-S观测**：

$$d_{NS} = N$$

**Vertical观测**：

$$d_{V} = U$$

#### 2.1.2 设计矩阵构建

对于局部窗口内的观测，建立观测方程：

$$L = B \cdot X + \Delta$$

其中：
- $L$：观测向量（$n \times 1$）
- $B$：设计矩阵（$n \times m$）
- $X$：未知参数向量（$m \times 1$）
- $\Delta$：观测误差（$n \times 1$）

**三维位移模型**（$fsmpara=2$，三维）：

未知参数包括：
- 3个位移参数：$[E, N, U]$
- 6个应变参数：$[\epsilon_{EE}, \epsilon_{EN}, \epsilon_{EU}, \epsilon_{NN}, \epsilon_{NU}, \epsilon_{UU}]$

设计矩阵为：

$$B_i = \begin{bmatrix} a_i & b_i & c_i & a_i \cdot de & a_i \cdot dn & b_i \cdot de & b_i \cdot dn & c_i \cdot de & c_i \cdot dn \end{bmatrix}$$

其中：
- $a_i, b_i, c_i$：第$i$种观测的几何系数
- $de, dn$：相对于中心像元的东向和北向坐标增量（米）

**三维位移模型**（$fsmpara=3$，三维，含高程）：

未知参数包括：
- 3个位移参数：$[E, N, U]$
- 9个应变参数：$[\epsilon_{EE}, \epsilon_{EN}, \epsilon_{EU}, \epsilon_{NE}, \epsilon_{NN}, \epsilon_{NU}, \epsilon_{UE}, \epsilon_{UN}, \epsilon_{UU}]$

设计矩阵为：

$$B_i = \begin{bmatrix} a_i & b_i & c_i & a_i \cdot de & a_i \cdot dn & a_i \cdot du & b_i \cdot de & b_i \cdot dn & b_i \cdot du & c_i \cdot de & c_i \cdot dn & c_i \cdot du \end{bmatrix}$$

其中$du$为相对于中心像元的高程增量。

### 2.2 加权最小二乘（WLS）

#### 2.2.1 基本公式

对于加权最小二乘问题：

$$\min_X (L - B \cdot X)^T \cdot P \cdot (L - B \cdot X)$$

其中$P$为权矩阵（对角矩阵）。

法方程为：

$$N \cdot X = U$$

其中：
- $N = B^T \cdot P \cdot B$（法方程系数矩阵）
- $U = B^T \cdot P \cdot L$（法方程常数项）

解为：

$$X = N^{-1} \cdot U = (B^T \cdot P \cdot B)^{-1} \cdot B^T \cdot P \cdot L$$

#### 2.2.2 方差估计

未知参数的方差-协方差矩阵为：

$$D_{XX} = \sigma_0^2 \cdot N^{-1}$$

其中$\sigma_0^2$为单位权方差，可由残差估计：

$$\hat{\sigma}_0^2 = \frac{V^T \cdot P \cdot V}{n - m}$$

其中$V = B \cdot X - L$为残差向量。

### 2.3 方差分量估计（VCE）

#### 2.3.1 基本原理

VCE用于估计不同类型观测的方差分量。假设有$k$类观测，每类观测的方差为$\sigma_i^2$（$i=1,2,...,k$）。

VCE的基本方程为：

$$S \cdot \hat{\theta} = W$$

其中：
- $S$：系数矩阵（$k \times k$）
- $\hat{\theta}$：方差分量向量（$k \times 1$），$\hat{\theta} = [\sigma_1^2, \sigma_2^2, ..., \sigma_k^2]^T$
- $W$：常数项向量（$k \times 1$）

#### 2.3.2 系数矩阵S的构建

系数矩阵$S$的元素计算如下：

**对角元素**：

$$S_{ii} = k_i - 2 \cdot \text{tr}(N^{-1} \cdot N_i) + \text{tr}((N^{-1} \cdot N_i)^2)$$

其中：
- $k_i$：第$i$类观测的有效样本数
- $N = \sum_{i=1}^{k} N_i$
- $N_i = B_i^T \cdot P_i \cdot B_i$

**非对角元素**：

$$S_{ij} = \text{tr}(N^{-1} \cdot N_i \cdot N^{-1} \cdot N_j)$$

#### 2.3.3 常数项W的计算

$$W_i = V_i^T \cdot P_i \cdot V_i$$

其中$V_i = B_i \cdot X - L_i$为第$i$类观测的残差向量。

#### 2.3.4 方差分量解算

$$\hat{\theta} = S^{-1} \cdot W$$

#### 2.3.5 权重更新

根据估计的方差分量更新各类观测的权重：

$$P_i^{new} = \frac{\hat{\sigma}_1^2}{\hat{\sigma}_i^2} \cdot P_i^{old}$$

其中$\hat{\sigma}_1^2$为参考类（第一类）的方差分量。

### 2.4 SM-VCE完整流程

SM-VCE的完整计算流程如下：

1. **初始化**：设置各类观测初始权重$P_i = I$（单位矩阵）
2. **第一次解算**：使用初始权重解算未知参数$X$
3. **方差分量估计**：根据残差估计各类观测的方差分量$\hat{\theta}$
4. **权重更新**：根据方差分量更新权重
5. **第二次解算**：使用更新后的权重重新解算未知参数
6. **输出结果**：输出最终的位移估计和方差估计

### 2.5 鲁棒VCE（可选）

鲁棒VCE用于抵抗粗差观测的影响，其权重调整函数为：

$$P_i^{robust} = \begin{cases} 0 & \text{if } |b_i| \geq t_2 \\ \frac{t_1 \cdot P_i}{|b_i|} \cdot \left(\frac{t_2 - |b_i|}{t_2 - t_1}\right)^2 & \text{if } t_1 \leq |b_i| < t_2 \\ P_i & \text{if } |b_i| < t_1 \end{cases}$$

其中：
- $b_i = \frac{|v_i|}{\sqrt{d_{tem} + \epsilon}}$
- $d_{tem} = \left(\frac{\text{median}(|v|)}{0.6745}\right)^2$
- $t_1, t_2$：阈值参数（默认$t_1=2, t_2=5$）

### 2.6 同质点筛选（SMAD）

SMAD（Strain Model Adaptive Neighbourhood）用于在局部窗口中剔除不满足平稳假设的样本。

#### 2.6.1 方向模板生成

生成32个方向模板，每个模板将窗口分为两个半区。方向角度为：

$$\theta_d = d \cdot \frac{360°}{32}, \quad d = 0, 1, ..., 31$$

#### 2.6.2 梯度计算

对于每个方向模板，计算梯度：

$$g_d = \sum_{i,j} smp_{ij} \cdot template_{ij}^d$$

其中$smp$为窗口内的观测值。

#### 2.6.3 主方向确定

找到使$|g_d|$最大的方向$d_{max}$，其垂直方向$d_{max}+90°$和$d_{max}-90°$为潜在的分界方向。

#### 2.6.4 线性拟合与筛选

1. 对窗口内观测进行线性拟合：$L = a + b \cdot x + c \cdot y$
2. 计算残差：$r = L_{obs} - L_{fit}$
3. 根据残差阈值筛选同质点：$|r| / \sigma \leq \sigma_{scale}$

其中$\sigma$为拟合残差的标准差，$\sigma_{scale}$为尺度因子（默认第一阶段$\sigma_{scale}=3$，第二阶段$\sigma_{scale}=6$）。

### 2.7 三维形变方差估计

对于解算得到的三维位移$X = [E, N, U]^T$，其方差-协方差矩阵为：

$$D_{ENU} = \hat{\sigma}_0^2 \cdot (B_{geo}^T \cdot P \cdot B_{geo})^{-1}$$

其中：
- $B_{geo}$：几何系数矩阵（$k \times 3$）
- $P$：权重对角矩阵
- $\hat{\sigma}_0^2$：单位权方差（取各类观测方差分量的平均值）

各分量的中误差为：

$$\sigma_E = \sqrt{|D_{ENU}[0,0]|}$$
$$\sigma_N = \sqrt{|D_{ENU}[1,1]|}$$
$$\sigma_U = \sqrt{|D_{ENU}[2,2]|}$$

---

## 3. 详细代码流程

### 3.1 主程序流程（smvce_main.py）

```
开始
  ↓
读取数据目录
  ↓
读取观测数据和几何信息
  ↓
是否运行WLS？ → 是 → 执行WLS三维反演
  ↓否
是否运行SM-VCE？ → 是 → 执行SM-VCE三维反演
  ↓否
绘制结果图
  ↓
导出结果
  ↓
结束
```

### 3.2 SM-VCE详细流程（smvce_solve3d.py）

```
开始
  ↓
初始化输出数组
  ↓
计算观测几何系数Bgeo
  ↓
检查是否可以使用CUDA快路径
  ↓
对于每个像元(i,j)：
  ├─ 检查mask
  ├─ 获取窗口观测数据_getL
  │   ├─ 提取窗口
  │   ├─ 断层约束（如有）
  │   └─ SMAD筛选（如启用）
  ├─ 构建观测向量和设计矩阵_getLBkP
  │   ├─ 计算坐标增量de, dn, du
  │   ├─ 构建设计矩阵B
  │   └─ 初始化权重P
  ├─ 几何秩和条件数检查
  ├─ 执行VCE：smvce_vce(L, B, P)
  │   ├─ 第一次解算_getsita
  │   ├─ 估计方差分量
  │   ├─ 更新权重_getP
  │   └─ 第二次解算_getsita
  ├─ 提取位移和应变参数
  └─ 估计三维形变方差
  ↓
输出结果enu, var_enu, para_sm, sita, SHPcount
  ↓
结束
```

### 3.3 VCE详细流程（smvce_vce.py）

```
开始
  ↓
初始化权重P
  ↓
第一次解算：_getsita(B, P, L)
  ├─ 构建法方程N = Σ(Bi^T·Pi·Bi)
  ├─ 构建常数项U = Σ(Bi^T·Pi·Li)
  ├─ 解算X = N^(-1)·U
  ├─ 计算残差Vi = Bi·X - Li
  ├─ 计算Wi = Vi^T·Pi·Vi
  ├─ 构建系数矩阵S
  └─ 解算方差分量θ = S^(-1)·W
  ↓
更新权重：Pi_new = (θ0/θi)·Pi_old
  ↓
第二次解算：_getsita(B, P_new, L)
  ↓
输出最终结果
  ↓
结束
```

---

## 4. 核心模块详解

### 4.1 lk_vec.py - 观测几何系数计算

**功能**：将各种观测类型的几何信息转换为对E/N/U的线性系数。

**输入**：
- `azi`：方位角数组（度）
- `inc`：入射角数组（度）
- `losazi`：观测类型标记（1-5）
- `leftright`：左右视标记

**输出**：
- `Bgeo`：几何系数数组（row × col × 3×data_num）

**公式**：
- LOS：$a = -\epsilon \sin\theta \sin\alpha$, $b = -\epsilon \sin\theta \cos\alpha$, $c = \cos\theta$
- AZI：$a = -\cos\alpha$, $b = \sin\alpha$, $c = 0$
- E-W：$a = 1$, $b = 0$, $c = 0$
- N-S：$a = 0$, $b = 1$, $c = 0$
- Vertical：$a = 0$, $b = 0$, $c = 1$

### 4.2 get_design_mat.py - 设计矩阵构建

**功能**：根据应变模型构建设计矩阵。

**函数**：
- `get_design_mat(Bgeo, de, dn, du)`：三维位移版本
- `get_design_mat_2D(Bgeo, de, dn, du)`：二维位移版本

**设计矩阵结构**：

三维模型（fsmpara=2）：
```
B = [a, b, c, a*de, a*dn, b*de, b*dn, c*de, c*dn]
```

三维模型（fsmpara=3）：
```
B = [a, b, c, a*de, a*dn, a*du, b*de, b*dn, b*du, c*de, c*dn, c*du]
```

### 4.3 smvce_vce.py - 方差分量估计

**功能**：实现VCE算法，估计各类观测的方差分量并更新权重。

**主要函数**：
- `smvce_vce(L, B, P)`：主函数
- `_getsita(B, P, L)`：计算方差分量
- `_getP(sita, v, P, neq0, t1, t2, vce_mode)`：更新权重
- `_robust_vce(P_vec, P_vec0, v, sita0, t1, t2)`：鲁棒VCE权重调整

### 4.4 smvce_smad.py - 同质点筛选

**功能**：基于SMAD方法筛选窗口内的同质点。

**主要函数**：
- `smvce_smad(smp, r0, c0, nom, iffault)`：主函数
- `_gettemps(siz, ndir)`：生成方向模板
- `_getminsum(smp, temps)`：基于梯度选择半区
- `_getmod33(smp)`：获取3×3块中值

**流程**：
1. 生成32个方向模板
2. 计算各方向梯度
3. 确定主方向和分界方向
4. 选择半区进行线性拟合
5. 根据残差阈值筛选同质点
6. 第二阶段精化筛选

### 4.5 wls3d.py - 加权最小二乘

**功能**：实现经典的WLS三维反演，作为基准参考。

**流程**：
1. 对每个像元，提取周围61×61窗口的几何系数
2. 构建法方程并解算
3. 输出位移和方差估计

### 4.6 smvce_solve3d.py - SM-VCE主求解器

**功能**：实现SM-VCE方法的核心求解逻辑。

**主要函数**：
- `smvce_solve3d(data_dict)`：CPU版本主函数
- `_smvce_solve3d_cuda_fixed_window(...)`：CUDA快路径
- `_getL(data, i, j, windowsize, coor, fault, flag_smad, flag_adpws)`：获取窗口观测
- `_getLBkP(...)`：构建观测向量、设计矩阵、权重矩阵
- `_getHomoPoints(...)`：断层同质点筛选

---

## 5. 输入输出约定

### 5.1 输入数据格式

**data_information文件**：每行5列，空格分隔
```
观测文件名 观测类型 入射角 方位角 左右视标记
```

**观测类型代码**：
- 1：LOS
- 2：AZI
- 3：E-W
- 4：N-S
- 5：Vertical

**可选文件**：
- `dem.tif`：数字高程模型
- `mask.tif`：处理掩膜（1=有效，0=无效）
- `fault.xy`：断层轨迹文件

### 5.2 输出结果

**输出目录**：`Result_SMVCE_时间戳/`

**输出文件**：
- `geotiff/*.tif`：GeoTIFF格式结果
- `grd/*.grd`：GMT GRD格式结果
- `README.md`：结果说明文件

**输出变量**：
- `enu_e`, `enu_n`, `enu_u`：三维位移
- `enu_std_e`, `enu_std_n`, `enu_std_u`：三维位移中误差
- `obs_std_*`：各类观测中误差
- `sita_*`：各类观测方差分量
- `shpcount_*`：各类观测有效样本数
- `para_sm_*`：应变模型参数

---

## 6. 使用方法

### 6.1 环境配置

```bash
# 激活conda环境
conda activate smvce_tiff

# 设置matplotlib缓存目录
export MPLCONFIGDIR=/tmp/mplconfig
```

### 6.2 运行程序

```bash
# 完整运行
MPLCONFIGDIR=/tmp/mplconfig conda run --no-capture-output -n smvce_tiff python smvce_main.py

# 或直接运行
python smvce_main.py
```

### 6.3 参数配置

在`smvce_main.py`顶部设置主要参数：

```python
# 求解控制
windowsize = 31          # 窗口大小
flag_if_2D = 0           # 0=三维，1=二维
fsmpara = 2              # 应变模型阶数（2或3）
flag_interWeight = 0     # 距离权重
flag_adpws = 0           # 自适应窗口
flag_smad = 0            # SMAD筛选
run_wls = 1              # 运行WLS
run_smvce = 1            # 运行SM-VCE

# CUDA参数
use_cuda = True          # 允许CUDA
cuda_rows_per_batch = 2  # CUDA批处理行数
```

---

## 7. 工程限制与注意事项

### 7.1 当前限制

1. **数据限制**：
   - 当前`SMVCE_DATA`中没有`dem.tif`，`fsmpara=3`时高程项不可靠
   - 当前`SMVCE_DATA`中没有`mask.tif`，所有像元都参与处理
   - 当前`SMVCE_DATA`中没有`fault.xy`，断层约束不生效

2. **环境限制**：
   - `shapely`未安装，断层几何判断受影响
   - CUDA不可用，默认使用CPU路径

3. **功能限制**：
   - `flag_interWeight`接口已保留但实现简化
   - 鲁棒VCE函数存在但默认未启用

### 7.2 性能考虑

1. **CPU路径**：逐像元循环，计算量大，耗时较长
2. **CUDA路径**：仅支持固定窗口+二维应变模型的默认配置
3. **窗口大小**：窗口越大，计算量越大，建议从31开始测试

### 7.3 结果验证

1. **WLS vs SM-VCE**：WLS结果作为基准参考，SM-VCE应给出更优的精度估计
2. **方差分量**：检查`sita`值是否合理，异常值可能表示VCE不收敛
3. **SHPcount**：检查有效样本数，过少可能表示窗口内数据不足

### 7.4 建议的改进方向

1. 安装`shapely`以启用断层约束
2. 配置CUDA环境以加速计算
3. 补充真实的`dem.tif`、`mask.tif`、`fault.xy`
4. 优化WLS模块的计算效率
5. 完善距离权重的实现

---

## 附录：参考文献

1. VCE理论基础：方差分量估计在大地测量中的应用
2. 应变模型：局部形变场的应变参数化方法
3. SMAD方法：基于应变模型的自适应邻域选择
4. InSAR观测几何：视线向、方位向形变的投影关系

---

## 8. 完整工作流程（2025定日 Mw6.8 地震案例）

本节描述从原始多源 InSAR 观测到最终三维形变场的完整处理流程。

### 8.1 输入数据组成

| 数据源 | 轨道 | 观测类型 | 文件 |
|--------|------|----------|------|
| Sentinel-1 | 升轨 | DInSAR LOS | `S1_Asc_DInSAR.disp.tif` |
| Sentinel-1 | 升轨 | POT LOS | `S1_Asc_POT_LOS.disp.tif` |
| LuTan-1 | 升轨 | DInSAR LOS | `LT1_Asc_DInSAR.disp.tif` |
| LuTan-1 | 升轨 | POT AZI | `LT1_Asc_POT_AZI.disp.tif` |
| LuTan-1 | 升轨 | MAI AZI | `LT1_Asc_MAI.disp.tif` |
| Sentinel-1 | 降轨 | DInSAR LOS | `S1_Des_DInSAR.disp.tif` |
| Sentinel-1 | 降轨 | POT LOS | `S1_Des_POT_LOS.disp.tif` |
| ALOS-2 | 降轨 | DInSAR LOS | `ALOS_Des_DInSAR.disp.tif` |

`data_information` 中按 **升降轨 → 卫星 → LOS/AZI** 顺序排列。

### 8.2 数据预处理

#### 8.2.1 LT1 DInSAR 常数偏移校正

LT1 DInSAR 与 S1 DInSAR 在远场存在系统性偏差（约 -0.33 m），原因是不同卫星的参考基准差异。

**校正方法**：

1. 选取远场区域（$\text{lon} < 87.38$ 或 $\text{lon} > 87.70$），排除形变区
2. 计算远场重叠区域 LT1 - S1 的中位差：

$$\Delta_{\text{offset}} = \text{median}\left(d_{\text{LT1}}^{\text{far}} - d_{\text{S1}}^{\text{far}}\right)$$

3. 施加常数校正：

$$d_{\text{LT1}}^{\text{corr}} = d_{\text{LT1}} - \Delta_{\text{offset}}$$

实测 $\Delta_{\text{offset}} \approx -0.328$ m，校正后远场残差均值 $< 0.07$ m。

#### 8.2.2 AZI 远场偏移对齐

POT AZI 与 MAI 观测在远场可能存在小偏移。以 MAI 为基准对齐 POT AZI：

$$d_{\text{POT}}^{\text{corr}} = d_{\text{POT}} - \text{median}\left(d_{\text{POT}}^{\text{far}} - d_{\text{MAI}}^{\text{far}}\right)$$

实测偏移约 -0.01 m，量级较小。

### 8.3 SM-VCE 反演

使用固定窗口 SM-VCE 方法（详见第 2 节）求解三维形变场 $[E, N, U]$。

**关键参数**：

| 参数 | 推荐值 | 说明 |
|------|--------|------|
| `windowsize` | 50 | 局部窗口尺寸（像素） |
| `fsmpara` | 2 | 应变模型阶数 |
| `flag_if_2D` | 0 | 三维求解 |
| `use_cuda` | True | CUDA 加速 |
| `cuda_rows_per_batch` | 自适应 | $\max(2, \lfloor 16 \times (40/ws)^2 \rfloor)$ |

CUDA batch 大小按窗口面积反比缩放，避免大窗口时显存溢出。

### 8.4 秩亏像素 2D 填充

在仅有 LOS 观测（无 AZI 覆盖）的区域，几何矩阵秩 < 3，无法求解三维。对这些像素使用升降轨 LOS 做 2D 分解（E-W + Vertical）：

$$\begin{bmatrix} d_{\text{asc}} \\ d_{\text{des}} \end{bmatrix} = \begin{bmatrix} c_E^{\text{asc}} & c_U^{\text{asc}} \\ c_E^{\text{des}} & c_U^{\text{des}} \end{bmatrix} \begin{bmatrix} E \\ U \end{bmatrix}$$

N-S 分量设为 0。

---

## 9. E-W 形变场后处理校正

SM-VCE 直接输出的 E-W 分量存在系统性偏差，主要表现为断层西侧出现不符合物理预期的正值（东向运动）。校正分两步进行。

### 9.1 Step 1: 远场平面去趋势（Deramp）

**目的**：消除 E-W 场中与经纬度相关的线性系统偏差。

**方法**：从远场像素（形变区外）拟合一阶平面趋势并移除。

1. 选取远场区域：$(\text{lon} < 87.38 \cup \text{lon} > 87.68) \cap (\text{lat} < 28.55 \cup \text{lat} > 29.00)$

2. 最小二乘拟合平面趋势：

$$E_{\text{ramp}}(\phi, \lambda) = a_0 + a_1 (\lambda - \bar{\lambda}) + a_2 (\phi - \bar{\phi})$$

其中 $\phi$ 为纬度，$\lambda$ 为经度，$\bar{\lambda}, \bar{\phi}$ 为全场中心。

3. 从全场移除趋势：

$$E_{\text{deramped}} = E_{\text{raw}} - E_{\text{ramp}}$$

实测趋势：$a_1 \approx 1.27$ m/°（经度方向）。

### 9.2 Step 2: 垂直分量泄漏校正（Vertical Leakage Correction）

**物理背景**：

在 LOS 到 ENU 的分解中，垂直分量（$U$）和东西分量（$E$）通过入射角强耦合。当升降轨 LOS 的几何多样性不足时，大幅垂直位移会"泄漏"到 E-W 估计中，产生与 $U$ 成正比的偏差。

**校正模型**：

$$E_{\text{corr}} = E_{\text{deramped}} + \lambda \cdot U$$

其中 $\lambda$ 为泄漏系数，$U$ 为垂直位移场（SM-VCE 输出，认为可靠）。

**$\lambda$ 的估算**：

定义断层两侧的近场区域：
- 西侧：$87.40° < \text{lon} < 87.52°$，$28.70° < \text{lat} < 28.90°$
- 东侧：$87.58° < \text{lon} < 87.68°$，$28.70° < \text{lat} < 28.90°$

**自动估算（反对称约束）**：

要求校正后西侧和东侧均值大小相等、符号相反：

$$\bar{E}_W^{\text{corr}} = -\bar{E}_E^{\text{corr}}$$

$$\bar{E}_W + \lambda \bar{U}_W = -(\bar{E}_E + \lambda \bar{U}_E)$$

解出：

$$\boxed{\lambda = -\frac{\bar{E}_W + \bar{E}_E}{\bar{U}_W + \bar{U}_E}}$$

实测自动估算 $\lambda \approx 1.46$，对应等效倾角 $\delta = \arctan(1/\lambda) \approx 34°$。

**注意**：自动估算的 $\lambda$ 通常过大（峰值校正量可达 2-3 m），原因是强制完美反对称的约束过强。**建议手动选择更保守的 $\lambda$**。

### 9.3 $\lambda$ 参数选择指南

| $\lambda$ | 西侧 E-W | 东侧 E-W | 峰值校正 | 适用场景 |
|-----------|----------|----------|---------|---------|
| 0（无校正） | +0.21 m | +0.46 m | 0 | 仅 deramp |
| 0.4 | ~0.00 m | +0.48 m | ~0.8 m | 保守：消除正值异常 |
| 0.7 | -0.16 m | +0.50 m | ~1.3 m | 折中：西侧出现负值 |
| 1.0 | -0.31 m | +0.52 m | ~1.9 m | 较强：西侧明显负值 |
| ~1.46（自动） | -0.55 m | +0.55 m | ~2.7 m | 强制反对称（过度校正风险） |

**推荐**：$\lambda \in [0.4, 0.7]$，结合正断层的预期 E-W 形变模式和已有研究结果选择。

### 9.4 校正的物理合理性

对于西倾正断层（定日地震，走向 ~167°，倾角 ~47°）：
- 上盘（西侧）：下沉 + 西移 → $U < 0$，$E < 0$
- 下盘（东侧）：微弱抬升 + 东移 → $U \approx 0$，$E > 0$

垂直位移场由升降轨 LOS 良好约束（$U_{\max} \approx -1.9$ m）。校正后 E-W 场应呈现**断层西侧负值（西移）、东侧正值（东移）**的模式。

---

## 10. 多剖面分析

SM-VCE 后处理自动生成 9 条东西向剖面（lat 28.55°N ~ 28.95°N，间隔 0.05°），保存至 `Result/Profile/` 目录。

每条剖面图包含 3 个子图：
- **(a)** SM-VCE 三维位移（E-W、N-S、Vertical）
- **(b)** 所有 LOS 原始观测
- **(c)** 所有 AZI 原始观测

`Profile_overview.png` 为总览图，展示 E-W 形变场及所有剖面线位置。

---

## 11. 命令行用法

```bash
# 基本运行（默认 ws=60）
conda run -n smvce_tiff python smvce_main.py

# 指定窗口大小
conda run -n smvce_tiff python smvce_main.py -w 50

# 添加标签
conda run -n smvce_tiff python smvce_main.py -w 50 --tag test

# 批量运行不同窗口大小
for ws in 40 50 60 80; do
    conda run -n smvce_tiff python -u smvce_main.py -w $ws
done
```

---

## 附录：参考文献

1. Liu, J. et al. SM-VCE method for 3D surface displacement estimation from multi-source InSAR observations.
2. Okada, Y. (1985). Surface deformation due to shear and tensile faults in a half-space. *BSSA*.
3. VCE 理论基础：方差分量估计在大地测量中的应用
4. SMAD 方法：基于应变模型的自适应邻域选择

---

*文档版本：v2.0*
*最后更新：2026年4月10日*
