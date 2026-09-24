"""Measure a version under an evaluation protocol.

    python -m vp_nrf2.evaluate --version v3

Writes ``versions/<version>/metrics.json`` and regenerates ``CARD.md``. The
protocol named in the version's manifest supplies the split, the fold sizes,
the seed set and the metric list. ``--protocol`` measures the same weights
under another registered protocol; ``metrics.json`` keys every run by
protocol id.

Each seed refits on its own training fold, so these are not the shipped
weights; see ``train``.

Every output is scored on the same folds, restricted to the compounds its own
endpoint labels. The first declared output carries the headline metrics.

A seed that yields an empty fold raises rather than being skipped.
"""

from __future__ import annotations

import argparse
import sys
from datetime import date

import numpy as np

from vp_core import fingerprints
from vp_nrf2 import are_ensemble
from vp_nrf2 import contract as nrf2_contract
from vp_nrf2 import data as nrf2_data
from vp_nrf2 import model as nrf2_model
from vp_nrf2.train import OUTPUT_LABELS, POOLED_OUTPUTS, pool_rows

__all__ = ["evaluate_version", "main"]


def _score_output(
    X, table, column, smiles, protocol, *, pool=None, features=None, epa_source=None
) -> tuple[dict, list[dict]]:
    """Per-seed metrics for one endpoint, on the folds its compounds fall in.

    When ``pool`` is given, its rows join the training fold, blocked by Murcko
    scaffold against the validation and test folds of that seed.
    """
    from vp_core import metrics as metrics_mod
    from vp_core.splits import murcko_scaffold

    scaffolds = np.array([murcko_scaffold(s) for s in smiles], dtype=object)
    rows = nrf2_data.labelled(table, column)
    keep = set(rows.tolist())
    y = table[column].to_numpy()

    per_seed: list[dict[str, float]] = []
    folds: list[dict[str, int]] = []
    for seed in protocol.seeds:
        train_idx, val_idx, test_idx = protocol.split_indices(smiles, seed)
        train = np.array([i for i in train_idx if i in keep])
        val = np.array([i for i in val_idx if i in keep])
        test = np.array([i for i in test_idx if i in keep])
        if not len(train) or not len(val) or not len(test):
            raise ValueError(
                f"seed {seed} leaves {column!r} with an empty fold — this endpoint "
                "cannot be measured under this protocol"
            )

        if epa_source is not None:
            blocked = set(scaffolds[np.concatenate([val, test])].tolist())
            fitted, n_extra = are_ensemble.fit_are(
                X, y.astype(int), train, val, blocked, epa_source, seed=seed
            )
            proba = are_ensemble.predict_are(fitted, X[test])
            scored = metrics_mod.binary_metrics(y[test].astype(int), proba)
            per_seed.append(scored)
            folds.append({"seed": int(seed), "train": len(train), "val": len(val),
                          "test": len(test), "external_train": n_extra})
            print(f"  {column} seed {seed}: test AUC {scored['auc_roc']:.4f} "
                  f"(n_test={len(test)}, n_external={n_extra})", file=sys.stderr)
            continue

        X_train, y_train = X[train], y[train].astype(int)
        n_extra = 0
        if pool is not None:
            blocked = set(scaffolds[val].tolist()) | set(scaffolds[test].tolist())
            blocked.discard("")
            extra_X, extra_y = pool_rows(pool, features, blocked=blocked)
            n_extra = len(extra_y)
            if n_extra:
                X_train = np.vstack([X_train, extra_X])
                y_train = np.concatenate([y_train, extra_y])

        fitted = nrf2_model.fit(X_train, y_train, X[val], y[val].astype(int), seed=seed)
        from vp_core import xgb

        proba = xgb.predict_proba(fitted, X[test])
        scored = metrics_mod.binary_metrics(y[test].astype(int), proba)
        per_seed.append(scored)
        folds.append(
            {
                "seed": int(seed),
                "train": len(train),
                "val": len(val),
                "test": len(test),
            }
        )
        print(
            f"  {column} seed {seed}: test AUC {scored['auc_roc']:.4f} "
            f"(n_test={len(test)})",
            file=sys.stderr,
        )

    return (
        {
            "n": {
                "total": len(rows),
                "train": folds[0]["train"],
                "val": folds[0]["val"],
                "test": folds[0]["test"],
            },
            "folds": folds,
            "test": metrics_mod.aggregate(per_seed, protocol.metrics),
        },
        folds,
    )


def evaluate_version(
    version: str,
    *,
    protocol_id: str | None = None,
    use_example: bool = False,
    write: bool = True,
) -> dict:
    """Run the version's protocol and return the metrics record."""
    import vp_nrf2
    from vp_core import card, dataset, metrics_store, protocols
    from vp_core import manifest as manifest_mod

    resolved = vp_nrf2.get(version)
    protocol_id = protocol_id or resolved.protocol_id
    protocol = protocols.get(protocol_id)
    labels = manifest_mod.dataset_labels(resolved.manifest)

    table = nrf2_data.example() if use_example else nrf2_data.load()
    smiles = table["smiles"].tolist()
    matrices: dict[str, np.ndarray] = {}

    def matrix(kind: str) -> np.ndarray:
        if kind not in matrices:
            matrices[kind] = fingerprints.featurize(smiles, kind)
        return matrices[kind]

    pooled = set(POOLED_OUTPUTS) & set(resolved.output_names)
    if pooled and not resolved.manifest.get("dataset", {}).get("auxiliary"):
        pooled = set()
    pool = None
    if pooled:
        pool = nrf2_data.example_pool() if use_example else nrf2_data.load_pool()
    epa_source = None
    if resolved.manifest["model"]["family"] == "xgboost-binary-ensemble":
        source_table = nrf2_data.example_epa() if use_example else nrf2_data.load_epa()
        epa_source = are_ensemble.source_features(
            are_ensemble.external_only(source_table, table)
        )

    scored = {
        output: _score_output(
            matrix(resolved.features_for(output)),
            table,
            OUTPUT_LABELS[output],
            smiles,
            protocol,
            pool=pool if output in pooled else None,
            features=resolved.features_for(output),
            epa_source=epa_source if output == "nrf2_are" else None,
        )[0]
        for output in resolved.output_names
    }
    primary = scored[nrf2_contract.PRIMARY]

    record = {
        "pathway": "nrf2",
        "version": resolved.name,
        "protocol_id": protocol_id,
        "dataset_sha256": dataset.dataset_hash(table, labels=labels),
        "dataset": "example fixture" if use_example else "full",
        "evaluated": date.today().isoformat(),
        "n": primary["n"],
        "folds": primary["folds"],
        "test": primary["test"],
        "additional_outputs": [
            {"name": name, **scored[name]}
            for name in resolved.output_names
            if name != nrf2_contract.PRIMARY
        ],
    }

    if write:
        if use_example:
            raise ValueError(
                "refusing to record fixture metrics as a released result — "
                "--example is for smoke-checking the harness only"
            )
        path = metrics_store.write(resolved.directory, record)
        card.write_card(resolved.directory)
        auc = record["test"]["auc_roc"]
        print(
            f"wrote {path}\n"
            f"  {protocol_id}: AUC {auc['mean']:.4f} +/- {auc['std']:.4f} "
            f"over {len(protocol.seeds)} seeds",
            file=sys.stderr,
        )
    return record


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m vp_nrf2.evaluate")
    parser.add_argument("--version", default=None, help="default: the newest version")
    parser.add_argument(
        "--protocol",
        default=None,
        help="default: the protocol the version's manifest declares",
    )
    parser.add_argument(
        "--example",
        action="store_true",
        help="run on the committed fixture without writing (smoke check only)",
    )
    args = parser.parse_args(argv)

    import vp_nrf2

    version = args.version or vp_nrf2.current_version()
    record = evaluate_version(
        version,
        protocol_id=args.protocol,
        use_example=args.example,
        write=not args.example,
    )
    if args.example:
        auc = record["test"]["auc_roc"]["mean"]
        print(f"fixture smoke AUC {auc if np.isfinite(auc) else float('nan'):.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
