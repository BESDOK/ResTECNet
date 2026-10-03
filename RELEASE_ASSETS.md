# Artefacts that are archived with each release but not stored in git

Git holds code and the paper tables only. The following are attached to the GitHub release and
archived on Zenodo (DOI in `CITATION.cff`), because of their size.

| Artefact | Source on the author's machine | Used by |
|---|---|---|
| Trained network weights (static, operational, uncertainty, Res-TECNet-A and the ablation/hyperparameter rows) | `*.mat` network files written by `runExperiments` (`resTECNet*.mat`, `ablation_*.mat`, `hyper_*.mat`) | `runExperiments` loads them instead of retraining |
| Processed station observations, one file per station and month | `stations/<CODE>_<yyyy>_<mm>.mat` | `evaluateStationVTEC`; the receiver-DCB harmonisation is reproduced by `harmonizeReceiverDCB` |
| Console log of the final run | `log_runAll.txt` | provenance |

Not archived (re-created by the code): the GIM and driver databases (`buildTECDatabase`,
`buildDriverTable`), the SH-AR forecast field `SHARforecasts.mat` (4.6 GB; `computeSHARForecasts`,
about 10 min), and the RINEX cache (`downloadRINEX`).
