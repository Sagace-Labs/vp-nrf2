# nrf2 v2

_Generated from `manifest.toml` and `metrics.json`. Do not edit._

Released 2026-09-19 · signature 1

**Why this version.** fits nrf2_cytotox on the panel's other viability counter-screens as well as its own

## Outputs

| column | dtype | range | meaning |
|---|---|---|---|
| `nrf2_are` | float32 | 0.0–1.0 | P(activates the antioxidant response element reporter in the Tox21 qHTS screen). The reporter registers indirect activators alongside direct electrophiles |
| `nrf2_cytotox` | float32 | 0.0–1.0 | P(reduces viability in the counter-screen over the same library). A compound scoring high on both readouts raised the reporter in a cell that was also dying |

Missing values: NaN when RDKit cannot parse the input SMILES

## Performance

Protocol `scaffold-balanced-5seed@1` — Bemis-Murcko scaffold split with scaffold groups permuted by seed and each group placed in the fold it overfills least, so a group larger than a fold's capacity settles in train instead of starving that fold. Same fold fractions, seeds and metrics as scaffold-shuffle-5seed@1; only the packing differs. Five seeds; report mean and standard deviation over the held-out test folds.

Evaluated 2026-09-19 on n_train=3496, n_val=417, n_test=627.

### `nrf2_are`

| metric | mean | std | per seed |
|---|---|---|---|
| auc_roc | 0.806 | 0.027 | 0.761, 0.845, 0.808, 0.814, 0.800 |
| auprc | 0.483 | 0.037 | 0.451, 0.531, 0.458, 0.449, 0.525 |
| mcc | 0.360 | 0.066 | 0.284, 0.427, 0.307, 0.332, 0.450 |
| brier | 0.127 | 0.044 | 0.212, 0.128, 0.097, 0.092, 0.106 |

### `nrf2_cytotox`

Measured on n_train=3400, n_val=406, n_test=595.

| metric | mean | std | per seed |
|---|---|---|---|
| auc_roc | 0.782 | 0.120 | 0.598, 0.908, 0.794, 0.700, 0.908 |
| auprc | 0.127 | 0.123 | 0.038, 0.355, 0.038, 0.045, 0.160 |
| mcc | 0.131 | 0.081 | 0.043, 0.220, 0.096, 0.063, 0.236 |
| brier | 0.169 | 0.070 | 0.243, 0.064, 0.210, 0.220, 0.107 |

> Comparable only with metrics carrying the same protocol id.

## Data

Tox21 antioxidant response element qHTS — PubChem BioAssay AID 743202 (qHTS assay for small molecule agonists of the antioxidant response element (ARE) signaling pathway) and AID 743203 (qHTS assay for small molecule agonists of the antioxidant response element (ARE) signaling pathway - cell viability counter screen), rows called Active or Inactive, one row per compound labelled by majority call across its assay records. Retrieved 2026-09-06, licensed public-domain, redistributed here.

`4540` compounds, positive rate `0.093`, table SHA-256 `9d966167f286084d…`

Regenerate and check for upstream drift with `python -m vp_nrf2.data fetch --verify`.

Tox21 viability counter-screens — PubChem BioAssay AID 743086 (qHTS assay to identify small molecule that activate the aryl hydrocarbon receptor (AhR) signaling pathway - cell viability counter screen), AID 720693 (qHTS assay to identify small molecule antagonists of the glucocorticoid receptor (GR) signaling pathway - cell viability counter screen), AID 743209 (qHTS assay for small molecule activators of the heat shock response signaling pathway - cell viability counter screen), AID 651633 (qHTS assay to identify small molecule agonists of the p53 signaling pathway - cell viability counter screen), AID 1346977 (qHTS assay to identify small molecule agonists of the pregnane X receptor (PXR) signaling pathway - cell viability counter screen), rows called Active or Inactive, one row per compound labelled by majority call across its assay records. Retrieved 2026-09-06, licensed public-domain, redistributed here.

`8177` compounds, table SHA-256 `5ea39f1bc029e196…`. Role: training-only.

## Model

xgboost-binary on `combo3` features. Shipped weights: one model per output, each on every compound its endpoint labels minus a 10% scaffold carve used for early stopping; nrf2_cytotox also trains on the further viability screens, scaffold-blocked against that carve

`weights.joblib` SHA-256 `05b232b551dc57c0…`

## Provenance

Environment: python 3.11.11, rdkit 2026.03.5, xgboost 3.2.0.

Reproducibility is to this dataset hash and this environment, not bit-exact: the sources are live endpoints and RDKit descriptor values move between releases.
