"""The NRF2 model: binary XGBoost on a Morgan/MACCS/descriptor composite.

Hyperparameters and feature choice; the fitting recipe comes from ``vp_core.xgb``.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from vp_core import fingerprints, xgb

__all__ = ["FEATURES", "FEATURES_BY_OUTPUT", "HYPERPARAMS", "features_for", "fit", "predict"]

FEATURES = "combo3"

#: Outputs fit with a featuriser other than ``FEATURES``.
FEATURES_BY_OUTPUT: dict[str, str] = {}


def features_for(output: str) -> str:
    """The featuriser ``output`` is fit with."""
    return FEATURES_BY_OUTPUT.get(output, FEATURES)


HYPERPARAMS: dict[str, Any] = {
    "n_estimators": 2000,
    "learning_rate": 0.05,
    "max_depth": 6,
    "min_child_weight": 1.0,
    "subsample": 0.8,
    "colsample_bytree": 0.8,
    "reg_lambda": 2.0,
    "reg_alpha": 0.05,
    "gamma": 0.1,
    "early_stopping_rounds": 40,
}


def fit(X_train, y_train, X_val, y_val, *, seed: int = 0) -> Any:
    """Fit one endpoint's binary model. The validation fold stops boosting."""
    return xgb.fit_binary(
        X_train,
        y_train,
        X_val,
        y_val,
        params=dict(HYPERPARAMS),
        seed=seed,
    )


def predict(model: Any, smiles: list[str], *, output: str | None = None) -> np.ndarray:
    """Positive-class probability for arbitrary SMILES."""
    kind = FEATURES if output is None else features_for(output)
    return xgb.predict_proba(model, fingerprints.featurize(smiles, kind))
