#!/bin/bash
# =============================================================================
# build_smvce_data.sh —— 由映射表组装一套新的 SM-VCE 输入数据集
#
# 用法:
#   ./build_smvce_data.sh [mapping.tsv] [out_dir]
# 默认:
#   mapping.tsv = mapping_smvce_data.tsv
#   out_dir     = SMVCE_DATA_new
#
# 行为:
#   1. 新建 out_dir (若已存在则先清空 .tif / data_information，保留以便重跑)
#   2. 从备份 SMVCE_DATA_backup_20260606 复用 几何(inc/azi) 及 dem/mask/fault
#   3. 逐行读 mapping:
#        source=REUSE     -> 从备份复制旧位移 tif
#        source=*.grd     -> 从 RAWGRD 用 gdalwarp 裁剪/转换到统一网格 (810x1170)
#        source=PENDING   -> 跳过并计入告警 (数据未就绪)
#   4. 重新生成 out_dir/data_information
#   5. 用 smvce_tiff 环境的 python 逐文件核验尺寸/网格一致性
#
# 统一目标网格 (与旧 SMVCE_DATA 完全一致):
#   region = 87.25 / 87.70 / 28.40 / 29.05   res = 2 arc-sec   ->  810 x 1170
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
MAPPING="${1:-$HERE/mapping_smvce_data.tsv}"
OUTDIR="${2:-$HERE/SMVCE_DATA_new}"
BACKUP="$HERE/SMVCE_DATA_backup_20260606"
RAWGRD="/media/jazz/SSD1/TJZ/#Post/2026.04.13 小论文final/Data/0605_Rawgrd"
PYBIN="/home/jazz/anaconda3/envs/smvce_tiff/bin/python"

# 统一目标网格
TE_W=87.25; TE_S=28.40; TE_E=87.70; TE_N=29.05
RES=0.000555555555556
RESAMP="near"   # 网格整数对齐 -> near = 精确裁剪，不插值、不侵蚀 NaN 边缘

echo "================ build_smvce_data ================"
echo "  mapping : $MAPPING"
echo "  out     : $OUTDIR"
echo "  backup  : $BACKUP"
echo "  rawgrd  : $RAWGRD"
echo "  grid    : -te $TE_W $TE_S $TE_E $TE_N  -tr $RES  ($RESAMP)"
echo "=================================================="

[ -f "$MAPPING" ] || { echo "ERR: 映射表不存在: $MAPPING"; exit 1; }
[ -d "$BACKUP" ]  || { echo "ERR: 备份目录不存在: $BACKUP"; exit 1; }

mkdir -p "$OUTDIR"
# 清理旧的位移 tif 与 data_information（保留几何复用文件以加速重跑亦可，这里全清求干净）
find "$OUTDIR" -maxdepth 1 -name '*.tif' -delete 2>/dev/null || true
rm -f "$OUTDIR/data_information"

# ---- 复用几何 / dem / mask / fault ----
echo "[1] 复用几何与辅助文件 ..."
shopt -s nullglob
for f in "$BACKUP"/*.inc.tif "$BACKUP"/*.azi.tif "$BACKUP"/dem.tif "$BACKUP"/mask.tif "$BACKUP"/fault.xy; do
  cp -a "$f" "$OUTDIR/"; echo "    + $(basename "$f")"
done
shopt -u nullglob

# ---- data_information 头 ----
DI="$OUTDIR/data_information"
{
  echo "#Displacement_Measurements	LOS_AZI_ENU	Inc_Angle	Azi_Angle	Left_or_Right_looking"
  echo "# 由 build_smvce_data.sh 依据 $(basename "$MAPPING") 生成"
} > "$DI"

n_reuse=0; n_conv=0; n_pending=0; pend_list=""

echo "[2] 逐组处理位移观测 ..."
# 读映射（跳过空行/注释）
while read -r tgt geom inc azi lr src _rest; do
  [ -z "${tgt:-}" ] && continue
  case "$tgt" in \#*) continue;; esac

  if [ "$src" = "REUSE" ]; then
    if [ -f "$BACKUP/$tgt" ]; then
      cp -a "$BACKUP/$tgt" "$OUTDIR/$tgt"
      printf "    [REUSE ] %-26s <- backup\n" "$tgt"
      n_reuse=$((n_reuse+1))
    else
      printf "    [MISS  ] %-26s 备份中无此文件!\n" "$tgt"; pend_list="$pend_list $tgt(missing-backup)"; n_pending=$((n_pending+1)); continue
    fi
  elif [ "$src" = "PENDING" ]; then
    printf "    [PEND  ] %-26s 等待你放入 0605_Rawgrd\n" "$tgt"
    pend_list="$pend_list $tgt"; n_pending=$((n_pending+1)); continue
  else
    if [ -f "$RAWGRD/$src" ]; then
      gdalwarp -overwrite -q -t_srs EPSG:4326 \
        -te $TE_W $TE_S $TE_E $TE_N -tr $RES $RES \
        -r $RESAMP -of GTiff -ot Float32 -dstnodata nan \
        "$RAWGRD/$src" "$OUTDIR/$tgt"
      printf "    [CONV  ] %-26s <- %s\n" "$tgt" "$src"
      n_conv=$((n_conv+1))
    else
      printf "    [MISS  ] %-26s 源文件不存在: 0605_Rawgrd/%s\n" "$tgt" "$src"
      pend_list="$pend_list $tgt(missing-src)"; n_pending=$((n_pending+1)); continue
    fi
  fi
  # 成功的行写入 data_information
  printf "%s\t%s\t%s\t%s\t%s\n" "$tgt" "$geom" "$inc" "$azi" "$lr" >> "$DI"
done < "$MAPPING"

echo "[3] 网格一致性核验 ..."
"$PYBIN" - "$OUTDIR" <<'PY'
import sys, os, glob
import rasterio, numpy as np
d=sys.argv[1]
# 以第一行 data_information 决定参考网格
di=os.path.join(d,'data_information')
disp=[]
with open(di) as f:
    for ln in f:
        ln=ln.strip()
        if not ln or ln.startswith('#'): continue
        disp.append(ln.split()[0])
ref=None; ok=True
print(f"  data_information 列出 {len(disp)} 组位移观测")
for nm in disp:
    p=os.path.join(d,nm)
    with rasterio.open(p) as s:
        a=s.read(1).astype(float)
        if s.nodata is not None: a[a==s.nodata]=np.nan
        key=(s.width,s.height,round(s.transform.c,7),round(s.transform.f,7),
             round(s.transform.a,12),round(s.transform.e,12))
        fin=np.isfinite(a)
        print(f"    {nm:28s} {s.width}x{s.height}  valid={100*fin.mean():4.1f}%  "
              f"range=[{np.nanmin(a):+.2f},{np.nanmax(a):+.2f}]")
        if ref is None: ref=key
        elif key!=ref: print(f"      !! 网格不一致: {key} != {ref}"); ok=False
print("  ==> 网格一致" if ok else "  ==> 警告: 存在网格不一致, SM-VCE 会堆叠错位!")
PY

echo "=================================================="
echo "  复用(REUSE)=$n_reuse  转换(CONV)=$n_conv  待处理(PENDING/MISS)=$n_pending"
[ -n "$pend_list" ] && echo "  待补:$pend_list"
echo "  输出: $OUTDIR"
echo "  下一步: 待全部就绪后, 备份并切换 ->"
echo "     mv SMVCE_DATA SMVCE_DATA_old_\$(date +%Y%m%d) && mv $OUTDIR SMVCE_DATA"
echo "=================================================="
