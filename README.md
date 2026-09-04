# SM-VCE group meeting materials

This repository contains the code and final slide deck used for the SM-VCE group meeting on 2026-09-04.

## Contents

- `matlab/`: Liu Jihong's MATLAB SM-VCE implementation and the group-meeting runner.
- `python/`: Python port and utilities used for the earthquake-case workflow.
- `presentation/`: the final 33-slide presentation with full-scene results.

## Data policy

Raw InSAR observations, ancillary rasters, fault files, generated grids, numerical results, figures, caches, and local environments are intentionally excluded. The repository therefore cannot reproduce the case-study figures without separately supplied data.

Place local input files in an ignored `SMVCE_DATA/` directory, or point the Python program to a local data directory. Do not commit those files.

## Python quick start

```powershell
cd python
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
python smvce_main.py --help
python smvce_main.py --data-dir D:\path\to\SMVCE_DATA --windowsize 41
```

The `--preliminary-gpu-no-fault-separation` option is only for fast full-scene preview. It disables cross-fault neighborhood removal and must not be treated as the strict reference solution.

## MATLAB quick start

Add `matlab/SMVCE_code/` and, when needed on newer MATLAB versions, `matlab/compatibility/` to the MATLAB path. Prepare a local `matlab/SMVCE_DATA/` folder and run `SMVCE_main.m`. The data directory is ignored by Git.

## Scope and attribution

The MATLAB implementation is the original version supplied by Liu Jihong. The Python version includes case-specific adaptations used in the presentation. No general-purpose parameter set is implied; window size, fault geometry, missing-data pattern, and stripe noise should be tuned for each earthquake case.

