# Reproducing the external measurements

The passive record is the LUMICKS C-Trap calibration deposited by Michael Bugiel on 14 March 2023. Use the file from version 7729824. The other record is the experimental trajectory in Example 1b from Argun and colleagues. The source notebook identifies the 50,000 samples as a measured harmonic-trap coordinate. The repository contains derived counts and scores. Download the raw files from the authors' archives. [C-Trap deposit](https://doi.org/10.5281/zenodo.7729824) · [C-Trap extraction example](https://lumicks-pylake.readthedocs.io/en/v1.6.1/tutorial/force_calibration/low_level_api.html) · [DeepCalib notebook and data](https://github.com/softmatterlab/DeepCalib/tree/bea01cc09413d9a19c5e23aad9e594fe879fea61)

From the project root in Windows PowerShell, use the P2 environment.

```powershell
New-Item -ItemType Directory -Force incoming | Out-Null
curl.exe -fL --retry 2 "https://zenodo.org/records/7729824/files/passive_calibration.h5" -o incoming/p2-passive.h5
curl.exe -fL --retry 2 "https://raw.githubusercontent.com/softmatterlab/DeepCalib/bea01cc09413d9a19c5e23aad9e594fe879fea61/Data_Example1b.mat" -o incoming/p2-harmonic.mat
.\.venv\Scripts\python.exe -m pip install -r empirical/requirements.txt
.\.venv\Scripts\python.exe -m unittest -v empirical.test_trajectory
.\.venv\Scripts\python.exe -m empirical.analyze --passive incoming/p2-passive.h5 --harmonic incoming/p2-harmonic.mat --out incoming/replay
```

On other platforms use the same Python module entry points. The reader checks the published HDF5 MD5 or the pinned MAT SHA-256 before analysis. Outputs include the SHA-256 of each original recording. The passive data use bead 1x and 1y in the ten-second acquisition window of calibration item 6, before that item's calibration was applied. Each channel is standardized on the training prefix, so absolute force calibration is unnecessary. Thresholds then remain fixed. All lags are measured in acquisition frames.

The two CSVs contain symmetric joint mass, floors, global Chapman–Kolmogorov prediction and directional discrepancies for every planned setting. Their JSON companions contain integer directional matrices at both lags, held-out sample sizes, all matched-visit statistics and environment metadata. The manifest records source and output hashes. Compare replay outputs with [the archived results](../results/trajectories/) to check independent reproduction. Two Windows runs with Python 3.12.10 and the pinned packages matched on all five artifacts.
