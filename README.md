# vp-nrf2

Predicts activation of the NRF2 / KEAP1 antioxidant response from a SMILES
string — a molecular initiating event for drug-induced liver injury, where an
electrophile modifies KEAP1's reactive cysteines, releases NRF2 and switches on
the antioxidant response element.

## Install

    pip install vp-nrf2

## Use

    from vp_nrf2 import predict
    predict(["CC(=O)Oc1ccccc1C(=O)O"])   # -> DataFrame[nrf2_are, nrf2_cytotox]

Returns one row per input and one column per declared output. `nrf2_are` is the
probability of raising the antioxidant response element reporter; `nrf2_cytotox`
is the probability of reducing viability in the counter-screen over the same
library. A compound scoring high on both raised the reporter in a cell that was
also dying, which the first column alone cannot say. The reporter registers
indirect activators alongside direct electrophiles. Unparseable SMILES come
back as NaN. Pin a version with `predict(smiles, version="v1")`; list what is
available with `versions()`.

## Current version

**v3**, signature 1, measured under protocol `scaffold-balanced-5seed@1`. The
full record — metrics per output and per seed, dataset hash, environment — is in
[`src/vp_nrf2/versions/v3/CARD.md`](src/vp_nrf2/versions/v3/CARD.md).

## Data

PubChem BioAssay AID 743202, the Tox21 qHTS screen for agonists of the
antioxidant response element signalling pathway, and AID 743203, its
cell-viability counter-screen over the same library. Both are reduced to one row
per compound labelled by the majority call across its assay records, retrieved
2026-09-06 and redistributed here as a United States government work in the
public domain. `nrf2_cytotox` also trains on five further Tox21 viability
counter-screens. `nrf2_are` also trains on EPA invitrodb v4.3 AEID 1110
calls, using archived EPA structure identities; both auxiliary tables are
training-only. Rebuild and check with `python -m vp_nrf2.data fetch --verify`
and `python -m vp_nrf2.data fetch-epa --verify`; see
[`data/README.md`](data/README.md) for the expected layout.

## Retrain

    python -m vp_nrf2.train --version v4 --reason "why this version exists" --supersedes v3 --recipe epa-assay-token
    python -m vp_nrf2.evaluate --version v4

`train` fits five calibrated reporter members and one viability model, with a
scaffold carve for early stopping and calibration. `evaluate` refits per seed
under the protocol. The folds were used to choose this recipe, so the reported
gain is exploratory. Reproducibility is to the recorded data and environment.

## Licence

Code is Apache-2.0 ([`LICENSE`](LICENSE)). The bundled dataset is in the public
domain ([`LICENSE-DATA`](LICENSE-DATA)).

## Cite

See [`CITATION.cff`](CITATION.cff).
