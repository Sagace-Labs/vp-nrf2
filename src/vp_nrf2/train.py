"""Build a new released version.

    python -m vp_nrf2.train --version v4 --reason "why this version exists" --recipe epa-assay-token

Writes ``versions/<version>/manifest.toml`` and ``weights.joblib``.

The shipped weights are fit on the whole dataset, with a scaffold carve
held out only for early stopping. The metrics in ``metrics.json`` come from a
different set of fits: the protocol refits per seed on its own train fold.
"""

from __future__ import annotations

import argparse
import platform
import shutil
import sys
from datetime import date
from pathlib import Path

import numpy as np

from vp_core import fingerprints
from vp_nrf2 import are_ensemble
from vp_nrf2 import contract as nrf2_contract
from vp_nrf2 import data as nrf2_data
from vp_nrf2 import model as nrf2_model
from vp_nrf2.target import CYTOTOX, TARGET, VIABILITY_POOL

__all__ = ["build_version", "main", "pool_rows"]

VERSIONS_DIR = Path(__file__).resolve().parent / "versions"

DEFAULT_PROTOCOL = "scaffold-balanced-5seed@1"

#: Which label column supplies each declared output.
OUTPUT_LABELS: dict[str, str] = {
    "nrf2_are": "label",
    "nrf2_cytotox": "cytotox",
}

#: Outputs that also train on the further viability screens in
#: ``nrf2_data.load_pool()``. Their label is unchanged; only the training rows
#: grow, and donor rows are blocked by Murcko scaffold against the fold that
#: judges the fit.
POOLED_OUTPUTS: tuple[str, ...] = ("nrf2_cytotox",)


def pool_rows(pool, features: str, *, blocked: set[str]):
    """Donor rows and their calls, one row per compound per screen that called it."""
    from vp_core.splits import murcko_scaffold

    pool_smiles = pool["smiles"].tolist()
    scaffolds = [murcko_scaffold(s) for s in pool_smiles]
    free = np.flatnonzero([s not in blocked for s in scaffolds])
    X = fingerprints.featurize(pool_smiles, features)

    blocks_X, blocks_y = [], []
    for column in nrf2_data.POOL_LABELS:
        values = pool[column].to_numpy(dtype="float64", na_value=np.nan)
        rows = free[np.isfinite(values[free])]
        if not len(rows):
            continue
        blocks_X.append(X[rows])
        blocks_y.append(values[rows].astype(int))
    if not blocks_X:
        return np.empty((0, X.shape[1]), dtype=np.float32), np.empty(0, dtype=int)
    return np.vstack(blocks_X), np.concatenate(blocks_y)


def _provenance() -> dict:
    """The environment that produced the version."""
    import rdkit
    import xgboost

    return {
        "python": platform.python_version(),
        "rdkit": rdkit.__version__,
        "xgboost": xgboost.__version__,
    }


def build_version(
    version: str,
    *,
    reason: str,
    protocol: str = DEFAULT_PROTOCOL,
    supersedes: str | None = None,
    seed: int = 0,
    recipe: str = "epa-assay-token",
) -> Path:
    """Fit the deployment models and write the version directory."""
    from vp_core import protocols
    from vp_core.splits import murcko_scaffold, scaffold_train_val

    protocols.get(protocol)  # fail early on an unknown protocol
    if recipe not in ("standard", "epa-assay-token"):
        raise ValueError(f"unknown NRF2 recipe: {recipe}")

    directory = VERSIONS_DIR / version
    if directory.exists():
        raise FileExistsError(
            f"{directory} already exists. Released versions are immutable — "
            "publish a new version instead of editing this one."
        )

    table = nrf2_data.load()
    smiles = table["smiles"].tolist()
    matrices: dict[str, np.ndarray] = {}

    def matrix(kind: str) -> np.ndarray:
        if kind not in matrices:
            matrices[kind] = fingerprints.featurize(smiles, kind)
        return matrices[kind]

    pool = nrf2_data.load_pool() if POOLED_OUTPUTS else None
    source = None
    if recipe == "epa-assay-token":
        source = are_ensemble.source_features(
            are_ensemble.external_only(nrf2_data.load_epa(), table)
        )

    fitted = {}
    for output in nrf2_contract.column_names():
        X = matrix(nrf2_model.features_for(output))
        column = OUTPUT_LABELS[output]
        rows = nrf2_data.labelled(table, column)
        y = table[column].to_numpy()[rows].astype(int)
        subset = [smiles[i] for i in rows]

        # Deployment fit: everything this endpoint labels, with a small scaffold
        # carve that stops boosting before it overfits.
        train_idx, val_idx = scaffold_train_val(subset, val_frac=0.10, seed=seed)
        if output == "nrf2_are" and source is not None:
            fitted[output], n_extra = are_ensemble.fit_are(
                X[rows], y, train_idx, val_idx,
                {murcko_scaffold(subset[i]) for i in val_idx}, source,
                seed=seed,
            )
            print(f"  {output}: {n_extra} EPA training rows, five members", file=sys.stderr)
            continue
        X_train, y_train = X[rows][train_idx], y[train_idx]

        if output in POOLED_OUTPUTS and pool is not None:
            extra_X, extra_y = pool_rows(
                pool,
                nrf2_model.features_for(output),
                blocked={murcko_scaffold(subset[i]) for i in val_idx} - {""},
            )
            X_train = np.vstack([X_train, extra_X])
            y_train = np.concatenate([y_train, extra_y])
            print(f"  {output}: {len(extra_y)} further viability rows", file=sys.stderr)

        fitted[output] = nrf2_model.fit(
            X_train,
            y_train,
            X[rows][val_idx],
            y[val_idx],
            seed=seed,
        )
        print(
            f"  {output}: {len(rows)} labelled compounds, "
            f"positive rate {y.mean():.3f}, "
            f"best iteration {getattr(fitted[output], 'best_iteration', None)}",
            file=sys.stderr,
        )

    directory.mkdir(parents=True)
    try:
        return _write_version(
            directory,
            fitted,
            table,
            version=version,
            reason=reason,
            protocol=protocol,
            supersedes=supersedes,
            recipe=recipe,
        )
    except Exception:
        shutil.rmtree(directory, ignore_errors=True)
        raise


def _write_version(
    directory: Path,
    fitted: dict,
    table,
    *,
    version: str,
    reason: str,
    protocol: str,
    supersedes: str | None,
    recipe: str,
) -> Path:
    import joblib

    import vp_core
    from vp_core import dataset, hashing, manifest

    weights_path = directory / "weights.joblib"
    joblib.dump(fitted, weights_path)

    labels = list(nrf2_data.LABELS)
    record = {
        "schema": manifest.SCHEMA_VERSION,
        "pathway": "nrf2",
        "version": version,
        "released": date.today().isoformat(),
        "supersedes": supersedes,
        "reason": reason,
        "signature": nrf2_contract.as_manifest_table(),
        "dataset": {
            "name": "Tox21 antioxidant response element qHTS",
            "source": (
                f"PubChem BioAssay AID {TARGET.pubchem_aid} ({TARGET.assay_name}) "
                f"and AID {CYTOTOX.pubchem_aid} ({CYTOTOX.assay_name}), rows called "
                f"Active or Inactive, one row per compound labelled by majority "
                f"call across its assay records"
            ),
            "url": (
                "https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/"
                f"{TARGET.pubchem_aid}/concise/CSV"
            ),
            "retrieved": "2026-09-06",
            "licence": "public-domain",
            "redistributable": True,
            "path": "data/nrf2_tox21.parquet",
            "labels": labels,
            "sha256": dataset.dataset_hash(table, labels=labels),
            "n_rows": len(table),
            "base_rate": round(float(table["label"].mean()), 6),
            "fetch": "python -m vp_nrf2.data fetch --verify",
        },
        "model": {
            "family": "xgboost-binary-ensemble" if recipe == "epa-assay-token" else "xgboost-binary",
            "features": nrf2_model.FEATURES,
            "fit": (
                ("five assay-token members for nrf2_are, trained with EPA AEID 1110 "
                 "rows excluded by primary identity and validation scaffold, "
                 "averaged and guarded-Platt calibrated on the validation carve; "
                 "one model for nrf2_cytotox, each using a 10% scaffold carve"
                 if recipe == "epa-assay-token" else
                 "one model per output, each on every compound its endpoint labels "
                 "minus a 10% scaffold carve used for early stopping")
                + (
                    f"; {', '.join(POOLED_OUTPUTS)} also trains on the further "
                    "viability screens, scaffold-blocked against that carve"
                    if POOLED_OUTPUTS
                    else ""
                )
            ),
            "weights": "weights.joblib",
            "sha256": hashing.sha256_file(weights_path),
        },
        "protocol": {
            "id": protocol,
            "provider": "vp-core",
            "core_version": vp_core.__version__,
        },
        "provenance": _provenance(),
    }

    if POOLED_OUTPUTS:
        pool_table = nrf2_data.load_pool()
        pool_labels = list(nrf2_data.POOL_LABELS)
        record["dataset"]["auxiliary"] = [
            {
                "name": "Tox21 viability counter-screens",
                "role": "training-only",
                "source": (
                    "PubChem BioAssay "
                    + ", ".join(
                        f"AID {e.pubchem_aid} ({e.assay_name})" for e in VIABILITY_POOL
                    )
                    + ", rows called Active or Inactive, one row per compound "
                    "labelled by majority call across its assay records"
                ),
                "retrieved": "2026-09-06",
                "licence": "public-domain",
                "redistributable": True,
                "path": "data/nrf2_viability_pool.parquet",
                "labels": pool_labels,
                "sha256": dataset.dataset_hash(pool_table, labels=pool_labels),
                "n_rows": len(pool_table),
                "fetch": "python -m vp_nrf2.data fetch --verify",
            }
        ]

    if recipe == "epa-assay-token":
        epa_table = nrf2_data.load_epa()
        record["dataset"]["auxiliary"].append({
            "name": "EPA invitrodb v4.3 AEID 1110 ARE activity calls",
            "role": "training-only",
            "source": (
                "EPA v4.3 PubChem export of 18 August 2026, Active/Inactive calls "
                "paired with archived 2018 EPA chemical identities; unresolved, "
                "inconclusive and primary connectivity identities excluded"
            ),
            "retrieved": "2026-09-23",
            "licence": "public-domain",
            "redistributable": True,
            "path": "data/nrf2_epa_aeid1110.parquet",
            "labels": ["label"],
            "sha256": dataset.dataset_hash(epa_table, labels=("label",)),
            "n_rows": len(epa_table),
            "fetch": "python -m vp_nrf2.data fetch-epa --verify",
        })

    if nrf2_model.FEATURES_BY_OUTPUT:
        record["model"]["features_by_output"] = [
            {"output": name, "kind": kind}
            for name, kind in nrf2_model.FEATURES_BY_OUTPUT.items()
        ]

    problems = manifest.validate(record, version_dir=directory)
    if problems:
        raise ValueError(f"refusing to write an invalid manifest: {'; '.join(problems)}")
    manifest.write(directory / "manifest.toml", record)

    print(
        f"wrote {directory}\n"
        f"  {len(table)} compounds, positive rate {table['label'].mean():.3f}\n"
        f"  next: python -m vp_nrf2.evaluate --version {version}",
        file=sys.stderr,
    )
    return directory


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m vp_nrf2.train")
    parser.add_argument("--version", required=True, help="new version name, e.g. v4")
    parser.add_argument("--reason", required=True, help="why this version exists")
    parser.add_argument("--protocol", default=DEFAULT_PROTOCOL)
    parser.add_argument("--supersedes", default=None)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--recipe", choices=("standard", "epa-assay-token"),
        default="epa-assay-token",
    )
    args = parser.parse_args(argv)

    build_version(
        args.version,
        reason=args.reason,
        protocol=args.protocol,
        supersedes=args.supersedes,
        seed=args.seed,
        recipe=args.recipe,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
