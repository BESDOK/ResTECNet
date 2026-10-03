# Res-TECNet

MATLAB code for *Storm-Time Prediction of Global Ionospheric TEC Maps Using a Residual CNN
with Spatiotemporal Encoding* (Advances in Space Research, AISR-D-26-01824, revised version).

## Quick start (final workflow of the revision)
```matlab
% 1. put earthdata_netrc (NASA Earthdata credentials, see below) in the working folder
% 2. verifyRevisionCoverage        % requirement-to-code check and synthetic smoke tests
% 3. runAll                        % station processing, then every table and figure (results/)
```
`runAll` calls `runStationProcessing` (RINEX download and VTEC processing, resumable) and
`runExperiments` (all tables and figures; finished networks and station months are cached).
The values reported in the paper are in `paper_tables/` (one CSV per table); a rerun
reproduces them up to the run-to-run variation of the network training.

## Versions and timestamps
The state of the code that belongs to a manuscript version is identified by a git tag and
its commit hash (`git rev-parse <tag>^{commit}`), by the GitHub release, and by the Zenodo
DOI of that release (see `CITATION.cff`). Large artefacts (trained weights, processed station
files) are listed in `RELEASE_ASSETS.md` and archived with the release. `CHECKSUMS.sha256`
lists the SHA-256 of every file of the release (`sha256sum -c CHECKSUMS.sha256`).

This code base implements **every experiment of the revised manuscript**, including the
data-download chain. Run `verifyRevisionCoverage` first (Part A checks the
requirement-to-function mapping, Part B runs download-free synthetic smoke tests),
then `runExperiments` (all tables/figures) or `runFullPaperSims` (May-2024 window only).

## Requirements
* MATLAB (developed and run with R2026a; R2024a or later should work), Deep Learning Toolbox (dlnetwork/trainnet/dlarray).
* Navigation Toolbox for the RINEX-based sections (`rinexread`, `rinexinfo`):
  `processStationVTEC`, `sppPositioning`.
* System tools: `curl` (built into Windows 10/11, macOS, Linux) and a NASA Earthdata
  account (free, https://urs.earthdata.nasa.gov). Put the credentials in a file named
  `earthdata_netrc` in the MATLAB working folder (one line:
  `machine urs.earthdata.nasa.gov login <user> password <pw>`); the downloaders pick
  it up automatically (or pass `'cddisNetrc', path`). Do not commit this file.
  CDDIS is the primary source for all GIMs (1998-present) and RINEX; AIUB is the
  fallback. `.Z` (Unix compress) files of years <= 2022 need `uncompress`/`gzip` or
  7-Zip (Windows: install to `C:\Program Files\7-Zip`, detected automatically);
  `CRX2RNX` (Hatanaka) is needed for RINEX observation files.
* GPU recommended (Res-TECNet 1999–2018 training ≈ hours; ConvLSTM/Transformer longer).

## Data sources (all downloaded automatically, cached locally)
| Product | Function | Archive |
|---|---|---|
| CODE final GIM (CODG / COD0OPSFIN) | `buildTECDatabase(...,'product','CODG')` | ftp.aiub.unibe.ch/CODE/yyyy |
| CODE 1-day predicted GIM (C1PG / COD0OPSPRD) | `'product','C1PG'` | ftp.aiub.unibe.ch/CODE/yyyy |
| UPC rapid GIM (UQRG, 15 min) | `'product','UQRG'` | cddis.nasa.gov (Earthdata) → chapman.upc.es fallback |
| Kp, F10.7 | `buildDriverTable` | kp.gfz-potsdam.de |
| Dst | `buildDriverTable` | wdc.kugi.kyoto-u.ac.jp (final→provisional→realtime) |
| NOAA day-1 Kp forecast (operational mode) | `buildDriverTable(...,'mode','operational')` | archived `yyyymmdd_3-day-forecast.txt` files in `index_cache/kp_forecast` |
| Satellite DCBs | `readCODEDCB` | ftp.aiub.unibe.ch/CODE/yyyy/P1P2yymm.DCB.Z |
| IGS RINEX obs/nav | `downloadRINEX` | cddis.nasa.gov → igs.bkg.bund.de fallback |
| Anatolian / eastern-Mediterranean set (AUT1, ISTA, MERS, NICO, BSHM, ARUC) and 16 global IGS stations (see `stationList.m`) | `downloadRINEX`, `processStationVTEC` | open IGS/EPN archives (CDDIS, BKG) — **no licensed TUSAGA-Aktif data are used**. ANKR, TUBI, NSSP, NYA1 and OHI3 are not available in the open archives and were replaced |

## Requirement → code mapping
| Review item | Implementation |
|---|---|
| R2.M1 independent GNSS validation, 4 regimes, other GIM | `stationList`, `downloadRINEX`, `readCODEDCB`, `processStationVTEC`, `evaluateStationVTEC` (Table 4), UQRG via `buildTECDatabase` + `evaluateTEC(...,'RefY',Yu)` |
| R2.M2 target-epoch drivers not available at 24-h lead | `buildDriverTable('mode','operational')`, `buildTECDataset('mode','operational')` (persisted F10.7, trailing mean, NOAA Kp forecast, last Dst); hindcast/operational rows in Tables 2, 3, 5, 6, 7 |
| R2.M3 proper baselines, same inputs | `persistenceBaseline`, `baselineC1PG`, `baselineSHAR` (SH degree 15, direct AR), `modelConvLSTM`/`initConvLSTM`, `modelTransformer`/`initTransformer`, `trainModelFunction` (same loss/optimizer/early stopping) |
| R2.M4 full CORS processing, several stations, gradients | `processStationVTEC` (L4 levelled to P4, MW + L4-rate slip tests, arcs ≥ 30 min, CODE satellite DCBs, daily receiver DCB, MSLM, IPP, 3×MAD screening), `regionalGradients` (Table 5) |
| R2.M5 slant errors + positioning | `evaluateSlantDelay`, `klobucharModel`, `sppPositioning` (Table 6) |
| R2.M6 uncertainty | `createResTECNet(...,'uncertainty',true)`, `lossTECNLL` (Eq. 7), `calibrationStats` (Sect. 4.9) |
| R2.m1 73rd column | `readIONEXTEC('dropDuplicateLon',true)` → W = 72 everywhere |
| R2.m2 2 h → 1 h interpolation | `interpolationInfluence` (decimation test + native-only training), `fillTECGaps` |
| R2.m3 hyperparameters | `hyperparameterSearch` (Table 1) |
| R2.m4 per-driver ablation | `driverAblation` (`buildTECDataset('zeroChannels',{...})`, zero padding, target block, operational) (Table 7) |
| PDF review 1 (forecast on Kayseri figure) | `visualizeMay2024('forecasts',...,'station',...)` |
| PDF review 2 (81-day mean needs future flux) | `buildDriverTable('meanType','trailing')` (default) |
| DOCX A2.41/A2.42 (colour-bar labels, new map figure, stations on error map) | `plotForecastMaps`, `plotErrorMap('stations',...)`, `visualizeMay2024` |

## Memory
The dataset struct does not store the input tensor X. `buildTECDataset` keeps the scaled
maps once (`ds.Ts`, single) plus per-sample driver rows and indices; `tecGetBatch(ds, idx)`
assembles X for a mini-batch on demand (used by `trainResTECNet`, `trainModelFunction`,
`evaluateTEC`). Memory per split is ~20 KB per hourly map (1999-2018 training: ~3.6 GB)
instead of ~0.57 MB per sample. Convert the database to single once after loading
(`TEC = single(TEC)`), as `runExperiments` does.

## Res-TECNet-A (adaptive base + learned correction)
The first full run (27 Sep) showed that the static network beats persistence and the other
learned baselines but not the adaptive SH-AR extrapolation or CODE's C1PG. The final model
therefore uses the SH-AR 24-h forecast (computed once for the whole record from past data
only, `computeSHARForecasts`, file `SHARforecasts.mat`) as input channel L+17 and as the
global-skip base (`createResTECNet(...,'baseChannel', L+17)`), so the network learns the
correction to an adaptive extrapolation instead of to persistence. `runExperiments`
(`do.adaptive = true`) warm-starts it from the static network; the ablation row
`persistBase` quantifies the gain of the adaptive base.

## Compute budget presets (`runExperiments`, `budget.preset`)
* `'fast'` (default, ~20 GPU-hours in total): training samples every 4 h (validation and
  test hourly), mini-batch 64, 25 epochs / patience 6 / LR halved every 8 epochs for the
  hindcast network and the plain CNN; the operational and the uncertainty networks are
  fine-tuned from the hindcast network (`transferWeights`, 10 epochs, LR 3e-4);
  ConvLSTM / Transformer 15 epochs; ablation (10 rows, all under one budget) and the
  hyperparameter table (two values per factor) on 2011-2018 with a 6-h stride and 10 epochs.
  Sect. 3.3 and Table 1 of the manuscript must state these settings.
* `'full'`: 1999-2018 with a 3-h stride, 40 epochs, all ablation rows and the full grid
  (several days of GPU time).

## Power-failure resilience (resume)
Everything long-running is checkpointed; after an interruption simply re-run the same
command and it continues:
* `buildTECDatabase`: each completed year is saved as `<outFile>_yYYYY.mat` and reused.
* `trainResTECNet`: with `'checkpointDir'` (set by `runExperiments` to `checkpoints/<tag>`)
  a checkpoint is written every epoch; a restart resumes from the newest one with the
  correct learning-rate stage. Delete the folder to retrain from scratch.
* `trainModelFunction` (ConvLSTM / Transformer): `'checkpointFile'` stores parameters,
  Adam moments and early-stopping state every epoch.
* `runExperiments`: every trained network / baseline / ablation row / hyperparameter row
  is saved to its own `.mat` and skipped when present; station processing is saved per
  station-month in `stations/`.

## Notes and caveats
* The code was written without access to a MATLAB session; `verifyRevisionCoverage`
  is the first thing to run — it exercises every function on synthetic data and reports
  PASS/FAIL. Fix any environment-specific issue it surfaces before the long runs.
* `processStationVTEC` uses the MSLM with **h = 450 km, α = 0.9782** as stated in the
  manuscript; CODE's own MSLM uses 506.7 km. Both are options (`'hShell'`, `'alpha'`);
  keep the manuscript and the code consistent.
* The NOAA archive URL template in `buildDriverTable/fetchKpForecast` is a placeholder;
  if the SWPC archive layout differs, place the archived `3-day-forecast.txt` files in
  `index_cache/kp_forecast/yyyymmdd_3-day-forecast.txt`. Epochs without an archived
  forecast use the empirical error model (σ = 0.9 Kp units) and are flagged in
  `drivers.KpFcSource`.
* `rinexread` returns epochs in GPS time; the ≤ 18 s offset is irrelevant for hourly
  matching. RINEX-2 observables are mapped by `rinexread`; C/A-only receivers need the
  P1C1 DCB (`readCODEDCB(...,'P1C1')` → `'dcbP1C1'`).
* The regional set uses only open-access stations because TUSAGA-Aktif data require a
  licence. MERS (Mersin/Erdemli, ~260 km from Kayseri) is the open station nearest to
  Kayseri and is used for the slant-delay, positioning and Fig.-5 overlays
  (`cfg.posStation`); gradients use NSSP−ISTA (east–west) and ANKR−NICO (north–south).
* IRI-2020 maps are not generated here (external Fortran/web service); sample them onto
  the grid at `dsTe.epochsY` and pass as `statsIRI` (see `runExperiments`, Sect. 5).
* The obsolete scripts `run_plotErrorMap.m` and `runSims.m` are superseded by
  `runFullPaperSims.m` / `runExperiments.m` and are no longer part of the release.
* All tables are written to `results/*.csv`; figures to `fig0.jpg` (error map, Fig. 6),
  `fig1–3.jpg` (May 2024, Figs. 3–5), `fig4.jpg` (forecast maps, Fig. 2).
