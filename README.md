# SM-VCE group meeting materials

This is the Git repository rooted at the local `2026.09.04 今日组会` working directory. Only the public release subset is tracked.

## Tracked contents

- `02_MATLAB原版/`: Liu Jihong's MATLAB SM-VCE implementation and the group-meeting runner.
- `03_Python版/`: the Python port and utilities used for the earthquake-case workflow.
- `05_PPT/SM-VCE_今日组会_整景结果版_20260904_v10.pptx`: the final 33-slide presentation.

## Excluded contents

Raw InSAR observations, ancillary rasters, fault files, generated grids, numerical results, figures, papers, the local presentation script, virtual environments, caches, and build artifacts are intentionally excluded by `.gitignore`.

The repository therefore cannot reproduce the case-study figures without separately supplied data. Place local inputs in an ignored `SMVCE_DATA/` directory and do not commit them.

## Python quick start

```powershell
cd '03_Python版\smvce_tjz'
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
python smvce_main.py --help
python smvce_main.py --data-dir D:\path\to\SMVCE_DATA --windowsize 41
```

The `--preliminary-gpu-no-fault-separation` option is only a fast full-scene preview. It disables cross-fault neighborhood removal and is not the strict reference solution.

## MATLAB quick start

Enter `02_MATLAB原版/SMVCE_MATLAB_原版`, add `SMVCE_code/` to the MATLAB path, prepare a local `SMVCE_DATA/` directory, and run `SMVCE_main.m`. The wrapper `02_MATLAB原版/run_groupmeeting_experiment.m` writes experiments outside the original source directory.

## Scope

The MATLAB implementation is the original version supplied by Liu Jihong. The Python version includes case-specific adaptations used in the presentation. Window size, fault geometry, missing-data pattern, and stripe noise must be reconsidered for every earthquake case.

