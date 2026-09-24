"""EPA ARE augmentation with an assay token and guarded calibration.

External calls add training rows only. Primary endpoint labels and inference
features remain the original Tox21 reporter and molecular structure.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from sklearn.linear_model import LogisticRegression

from vp_core import fingerprints, xgb
from vp_core.splits import murcko_scaffold
from vp_nrf2 import model

MEMBERS = 5


@dataclass
class AreEnsemble:
    members: list
    slope: float | None
    intercept: float | None
    prior_weight: float


def external_only(source, primary):
    """Exclude primary connectivity identities and conflicting source stereoisomers."""
    primary_keys = set(primary.inchikey.str.split("-").str[0])
    out = source.loc[~source.inchikey.str.split("-").str[0].isin(primary_keys)].copy()
    out["block"] = out.inchikey.str.split("-").str[0]
    consistent = out.groupby("block").label.nunique()
    out = out.loc[out.block.isin(consistent[consistent == 1].index)]
    return out.drop_duplicates("block").sort_values("inchikey").reset_index(drop=True)


def source_features(source):
    """The standardized source's features and scaffold keys, computed once."""
    return (
        fingerprints.featurize(source.smiles.tolist(), model.FEATURES),
        source.label.to_numpy(dtype=int),
        np.asarray([murcko_scaffold(s) for s in source.smiles]),
    )


def _logit(values):
    p = np.clip(np.asarray(values, dtype=float), 1e-7, 1 - 1e-7)
    return np.log(p / (1 - p))


def fit_are(X, y, train, val, blocked, source, *, seed: int) -> tuple[AreEnsemble, int]:
    """Fit matched members and calibrate only on the validation fold."""
    source_X, source_y, source_scaffolds = source
    allowed = np.flatnonzero([s not in blocked for s in source_scaffolds])
    Xt = np.concatenate(
        [np.column_stack([X[train], np.zeros(len(train))]),
         np.column_stack([source_X[allowed], np.ones(len(allowed))])]
    )
    yt = np.concatenate([y[train], source_y[allowed]])
    Xv = np.column_stack([X[val], np.zeros(len(val))])
    fitted = []
    validation = []
    for member in range(MEMBERS):
        booster = xgb.fit_binary(
            Xt, yt, Xv, y[val], params=dict(model.HYPERPARAMS),
            seed=seed + 100 * member,
        )
        fitted.append(booster)
        validation.append(xgb.predict_proba(booster, Xv))
    mean_validation = np.mean(validation, axis=0)
    prior_weight = float((yt == 0).sum() / (yt == 1).sum())
    calibrator = LogisticRegression(C=1e6, max_iter=2000)
    calibrator.fit(_logit(mean_validation)[:, None], y[val])
    slope = float(calibrator.coef_[0, 0])
    if slope <= 0:
        return AreEnsemble(fitted, None, None, prior_weight), len(allowed)
    return AreEnsemble(
        fitted, slope, float(calibrator.intercept_[0]), prior_weight
    ), len(allowed)


def predict_are(fitted: AreEnsemble, X) -> np.ndarray:
    """Predict the primary endpoint with the source token fixed to zero."""
    primary = np.column_stack([X, np.zeros(len(X))])
    raw = np.mean([xgb.predict_proba(member, primary) for member in fitted.members], axis=0)
    if fitted.slope is None:
        return raw / (raw + fitted.prior_weight * (1 - raw))
    z = fitted.slope * _logit(raw) + fitted.intercept
    return 1 / (1 + np.exp(-z))
