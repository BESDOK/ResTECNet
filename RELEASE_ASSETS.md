# Artefacts that are archived with each release but not stored in git

Git holds code and the paper tables only. The files below are uploaded directly to Zenodo (DOI in `CITATION.cff`),
because of their size; Zenodo's GitHub integration archives only the source code of a tag, not files attached to a
release. File names are those written by the code
(`runExperiments`, `driverAblation`, `hyperparameterSearch`, `processStationYears`); they sit in the MATLAB
working folder unless a sub-folder is given.

## A. Trained network weights (MAT files; 26 files = 9 main + 9 ablation + 8 hyperparameter)

| File | Variable | Model in the paper |
|---|---|---|
| `resTECNet_hindcast.mat` | `net` | Res-TECNet, persistence base, hindcast (Tables 2, 3) |
| `resTECNet_operational.mat` | `netOp` | Res-TECNet, persistence base, operational |
| `resTECNet_uncertainty.mat` | `netU` | persistence base, heteroscedastic head |
| `resTECNetA_hindcast.mat` | `netA` | **Res-TECNet-A, hindcast (final model)** |
| `resTECNetA_operational.mat` | `netOpA` | **Res-TECNet-A, operational** |
| `resTECNetA_uncertainty.mat` | `netUA` | Res-TECNet-A, uncertainty head (Sect. 4.9) |
| `baseline_plain.mat` | `netPlain` | plain CNN (Table 2) |
| `baseline_convlstm.mat` | `pCL`, `cfgCL`, `infoCL` | ConvLSTM (Table 2) |
| `baseline_transformer.mat` | `pTR`, `cfgTR`, `infoTR` | Transformer (Table 2) |

Ablation rows (Table 7; fine-tuned from `resTECNetA_hindcast.mat`; tag `_A2`; `noisyTarget` is inference-only and has no file):

`ablation_A2_full.mat`, `ablation_A2_persistBase.mat`, `ablation_A2_noGlobalSkip.mat`, `ablation_A2_zeroPad.mat`,
`ablation_A2_noDrivers.mat`, `ablation_A2_noTargetBlock.mat`, `ablation_A2_noKp.mat`, `ablation_A2_noDst.mat`,
`ablation_A2_operational.mat`  (9 files)

Hyperparameter rows (Table 1; tag `_A2`):

`hyper_A2_lambda_0.mat`, `hyper_A2_lambda_0.05.mat`, `hyper_A2_lambda_0.1.mat`, `hyper_A2_lambda_0.2.mat`,
`hyper_A2_lambda_0.5.mat`, `hyper_A2_B_4.mat`, `hyper_A2_B_12.mat`, `hyper_A2_B_16.mat`  (8 files)

Do **not** include `ablation_A2_plainCNN.mat` (left over from an earlier budget; the plain CNN of Table 2 is `baseline_plain.mat`), nor superseded files from earlier protocols: any `ablation_A_*.mat`, `ablation_*` without the `_A2` tag,
any `hyper_P2_*.mat`, and anything under `checkpoints/`.
The network of the native-only experiment (Sect. 4.8.1) is not stored separately: `results/interpolation_influence.mat`
holds only its test statistics (and the 0.7 GB forecast array), and the text summary
`results/interpolation_influence_summary.txt` is archived instead.

## B. Processed station observations (510 monthly files + 1)

`stations/<CODE>_<yyyy>_<mm>.mat` (variable `obsAll`, one row per ionospheric pierce-point observation, before the
receiver-DCB harmonisation) for 2023-01 ... 2024-12, and `station_MERS.mat` (positioning station, `cfg.posStation`).

* 16 IGS stations: BOGT, IISC, NKLG, GUAM, ONSA, ALGO, MIZU, STR1, KIR0, THU2, MAC1, DAV1, KOKB, ASCG, DGAR, THTI
* 6 regional stations: AUT1, ISTA, MERS, NICO, BSHM, ARUC

Months that do not exist because no usable data were available in the archives (18 station-months; listed in
`stations_missing_months.csv`, also archived):

* ALGO: 2023-05, 2023-06, 2023-07, 2023-08
* MAC1: 2024-05
* KOKB: 2024-10, 2024-11
* AUT1: 2023-01, 2023-02, 2023-03, 2023-04, 2024-02, 2024-03
* ISTA: 2024-06
* BSHM: 2024-06, 2024-07, 2024-08, 2024-09

Expected count: 22 stations x 24 months - 18 = 510 files in `stations/`.

## C. Provenance and small inputs

`log_runAll.txt` (console log of the final run), `stations_missing_months.csv`, `drivers_hindcast.mat` (2.5 MB) and
`drivers_operational.mat` (5.7 MB) (the driver tables built by `buildDriverTable`), and
`results/interpolation_influence_summary.txt`.

## Not archived (re-created by the code)

GIM and driver databases (`TECdatabase.mat`, `C1PGdatabase.mat`, `UQRGdatabase.mat`; `buildTECDatabase`), the SH-AR forecast
field `SHARforecasts.mat` (4.6 GB; `computeSHARForecasts`, about 10 min), the RINEX, IONEX and DCB caches, and `checkpoints/`.

## Verifying the list on your machine (PowerShell, in the MATLAB working folder)

```powershell
Get-ChildItem resTECNet*.mat, baseline_*.mat, ablation_A2_*.mat, hyper_A2_*.mat | Select-Object Name, @{n='MB';e={[math]::Round($_.Length/1MB,1)}}
(Get-ChildItem stations\*.mat).Count            # expected: 510
Get-ChildItem station_MERS.mat, log_runAll.txt, stations_missing_months.csv, drivers_*.mat
```

## Checked against the author's disk (2 Oct 2026)
27 network files matched the patterns: the 26 listed above plus `ablation_A2_plainCNN.mat`, which is excluded.
Sizes: weights 1.4-17.6 MB each (about 70 MB in total; `baseline_transformer.mat` 17.6 MB), `stations/` 510 files,
`station_MERS.mat` 185 MB, `log_runAll.txt` 1.8 MB (it contains the console output of several runs, including
the interrupted ones). Zenodo accepts up to 50 GB per record; GitHub release assets are limited to 2 GB per file.

## Packaging (PowerShell; `Compress-Archive` fails above 2 GB, Windows `tar` does not)
```powershell
$w = (Get-ChildItem resTECNet*.mat, baseline_*.mat, ablation_A2_*.mat, hyper_A2_*.mat | Where-Object Name -ne 'ablation_A2_plainCNN.mat').Name
tar -a -cf ResTECNet_weights.zip $w
foreach ($y in 2023,2024) { tar -a -cf ("ResTECNet_stations_$y.zip") -C stations (Get-ChildItem stations\*_${y}_*.mat).Name }
tar -a -cf ResTECNet_provenance.zip station_MERS.mat log_runAll.txt stations_missing_months.csv drivers_hindcast.mat drivers_operational.mat results\interpolation_influence_summary.txt
(Get-ChildItem stations -Recurse | Measure-Object Length -Sum).Sum / 1GB     # total size of the station files (GB)
```
