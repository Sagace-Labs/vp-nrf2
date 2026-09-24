"""The NRF2 dataset: obtain, standardise, verify.

The shipped table is the parsed result — one row per compound, keyed on the
standardised InChIKey — and that is what the manifest hash covers. ``fetch``
rebuilds it from PubChem and compares hashes, so a differing hash means the
upstream assay record changed.

Two endpoints share the table. ``label`` is the antioxidant-response call and
defines which compounds the table holds. ``cytotox`` is the viability call from
the counter-screen run on the same library, and is null for a compound that
screen did not call.

Pipeline, per endpoint: pull the concise BioAssay table, keep the rows the assay
called Active or Inactive, resolve each compound identifier to a SMILES string
through the compound property endpoint, standardise, then collapse to one row
per compound by majority vote across its assay records.

Run as a command:

    python -m vp_nrf2.data fetch --verify
    python -m vp_nrf2.data verify
"""

from __future__ import annotations

import argparse
import io
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from vp_nrf2.target import CYTOTOX, TARGET, VIABILITY_POOL

__all__ = [
    "DATA_DIR",
    "EXAMPLE_PATH",
    "LABELS",
    "POOL_EXAMPLE_PATH",
    "POOL_LABELS",
    "POOL_PATH",
    "TABLE_PATH",
    "build_example",
    "example",
    "example_epa",
    "example_pool",
    "fetch",
    "labelled",
    "load",
    "load_epa",
    "load_pool",
    "verify",
]

_PACKAGE_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = _PACKAGE_ROOT / "data"
TABLE_PATH = DATA_DIR / "nrf2_tox21.parquet"
EXAMPLE_PATH = DATA_DIR / "example" / "nrf2_example.parquet"

#: The label columns the table carries, primary endpoint first.
LABELS: tuple[str, ...] = ("label", "cytotox")

POOL_PATH = DATA_DIR / "nrf2_viability_pool.parquet"
POOL_EXAMPLE_PATH = DATA_DIR / "example" / "nrf2_viability_pool_example.parquet"
EPA_PATH = DATA_DIR / "nrf2_epa_aeid1110.parquet"
EPA_EXAMPLE_PATH = DATA_DIR / "example" / "nrf2_epa_example.parquet"

#: One call column per further viability screen, in ``VIABILITY_POOL`` order.
POOL_LABELS: tuple[str, ...] = tuple(f"y_{e.name.lower()}" for e in VIABILITY_POOL)

PUG_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
_MIN_INTERVAL = 0.25  # PubChem asks for no more than five requests a second
_CID_BATCH = 200  # the property endpoint carries its identifiers in the URL

# The two calls that carry a label. A qHTS curve the assay could not call
# either way is not a negative, so "Inconclusive" is dropped.
_ACTIVE = "Active"
_INACTIVE = "Inactive"


def _missing(path: Path) -> str:
    return (
        f"{path} is not present. The dataset is not included in the wheel; run "
        "`python -m vp_nrf2.data fetch` to rebuild it from PubChem."
    )


def load() -> pd.DataFrame:
    """The full standardised dataset. Raises with guidance when absent."""
    from vp_core import dataset

    if not TABLE_PATH.exists():
        raise FileNotFoundError(_missing(TABLE_PATH))
    return dataset.read_table(TABLE_PATH, labels=LABELS)


def example() -> pd.DataFrame:
    """The committed test fixture — small, stratified, always available."""
    from vp_core import dataset

    if not EXAMPLE_PATH.exists():
        raise FileNotFoundError(_missing(EXAMPLE_PATH))
    return dataset.read_table(EXAMPLE_PATH, labels=LABELS)


def load_pool() -> pd.DataFrame:
    """The further viability screens. Raises with guidance when absent."""
    from vp_core import dataset

    if not POOL_PATH.exists():
        raise FileNotFoundError(_missing(POOL_PATH))
    return dataset.read_table(POOL_PATH, labels=POOL_LABELS)


def load_epa() -> pd.DataFrame:
    """EPA v4.3 ARE calls used only to augment reporter training."""
    from vp_core import dataset

    if not EPA_PATH.exists():
        raise FileNotFoundError(
            f"{EPA_PATH} is not present; run `python -m vp_nrf2.data fetch-epa`."
        )
    return dataset.read_table(EPA_PATH, labels=("label",))


def example_epa() -> pd.DataFrame:
    """A small stratified EPA fixture for the evaluation smoke run."""
    from vp_core import dataset

    if not EPA_EXAMPLE_PATH.exists():
        raise FileNotFoundError(
            f"{EPA_EXAMPLE_PATH} is not present; run `python -m vp_nrf2.data build-example`."
        )
    return dataset.read_table(EPA_EXAMPLE_PATH, labels=("label",))


def example_pool() -> pd.DataFrame:
    """The committed fixture for the further viability screens."""
    from vp_core import dataset

    if not POOL_EXAMPLE_PATH.exists():
        raise FileNotFoundError(_missing(POOL_EXAMPLE_PATH))
    return dataset.read_table(POOL_EXAMPLE_PATH, labels=POOL_LABELS)


def labelled(table: pd.DataFrame, column: str) -> np.ndarray:
    """Row positions where ``column`` carries a call."""
    return np.flatnonzero(table[column].notna().to_numpy())


# ---------------------------------------------------------------------------
# Rebuild from source
# ---------------------------------------------------------------------------


def _get(url: str, *, timeout: int, session=None):
    """GET with backoff. PubChem answers 503 when a rebuild asks too fast."""
    import requests

    getter = session.get if session is not None else requests.get
    delay = 1.0
    for attempt in range(6):
        response = getter(url, headers={"Accept": "text/csv"}, timeout=timeout)
        if response.status_code not in (429, 503) or attempt == 5:
            response.raise_for_status()
            return response
        time.sleep(delay)
        delay *= 2
    raise RuntimeError("unreachable")


def _assay_records(aid: int) -> pd.DataFrame:
    """The BioAssay table for one AID."""
    response = _get(f"{PUG_BASE}/assay/aid/{aid}/concise/CSV", timeout=300)
    raw = pd.read_csv(io.BytesIO(response.content), dtype=str, low_memory=False)

    for column in ("CID", "Activity Outcome"):
        if column not in raw.columns:
            raise RuntimeError(
                f"AID {aid} returned no {column!r} column; the concise format "
                f"changed upstream. Columns: {raw.columns.tolist()}"
            )

    df = pd.DataFrame(
        {
            "cid": pd.to_numeric(raw["CID"], errors="coerce"),
            "outcome": raw["Activity Outcome"].astype(str).str.strip(),
            "potency_um": pd.to_numeric(raw.get("Activity Value [uM]"), errors="coerce"),
        }
    )
    df = df[df["outcome"].isin([_ACTIVE, _INACTIVE])].dropna(subset=["cid"])
    df["cid"] = df["cid"].astype(int)
    print(
        f"  AID {aid}: {len(df)} labelled records "
        f"({int((df['outcome'] == _ACTIVE).sum())} active)",
        file=sys.stderr,
    )
    return df.reset_index(drop=True)


def _smiles_for(cids: list[int]) -> pd.DataFrame:
    """``(cid, smiles_raw)`` from the compound property endpoint, in batches."""
    import requests

    session = requests.Session()
    frames: list[pd.DataFrame] = []
    unique = sorted(set(cids))
    for start in range(0, len(unique), _CID_BATCH):
        batch = unique[start : start + _CID_BATCH]
        ids = ",".join(str(c) for c in batch)
        url = f"{PUG_BASE}/compound/cid/{ids}/property/SMILES,ConnectivitySMILES/CSV"
        response = _get(url, timeout=120, session=session)
        page = pd.read_csv(io.BytesIO(response.content))
        primary = "SMILES" if "SMILES" in page.columns else "IsomericSMILES"
        fallback = (
            "ConnectivitySMILES"
            if "ConnectivitySMILES" in page.columns
            else "CanonicalSMILES"
        )
        frames.append(
            pd.DataFrame(
                {
                    "cid": pd.to_numeric(page["CID"], errors="coerce"),
                    "smiles_raw": page[primary].fillna(page.get(fallback)),
                }
            )
        )
        print(
            f"  resolved {min(start + _CID_BATCH, len(unique))} / {len(unique)}",
            file=sys.stderr,
        )
        time.sleep(_MIN_INTERVAL)

    if not frames:
        raise RuntimeError("PubChem returned no compound structures at all.")
    out = pd.concat(frames, ignore_index=True).dropna(subset=["cid", "smiles_raw"])
    out["cid"] = out["cid"].astype(int)
    return out.drop_duplicates(subset=["cid"]).reset_index(drop=True)


def _standardise(records: pd.DataFrame, structures: pd.DataFrame) -> pd.DataFrame:
    from rdkit import Chem, RDLogger

    from vp_core.standardise import standardise_many

    RDLogger.DisableLog("rdApp.*")

    df = records.merge(structures, on="cid", how="inner")
    smiles, keys = standardise_many(df["smiles_raw"].astype(str).tolist())
    df["smiles"] = smiles
    df["inchikey"] = keys
    df = df.dropna(subset=["smiles", "inchikey"])

    # A featuriser turns a SMILES it cannot read into an all-zero row rather
    # than an error, so the round trip is required here.
    df = df[[Chem.MolFromSmiles(s) is not None for s in df["smiles"]]]

    df["active"] = (df["outcome"] == _ACTIVE).astype(int)
    return df


def _to_compounds(records: pd.DataFrame, structures: pd.DataFrame) -> pd.DataFrame:
    """Collapse reporter records to one row per standardised compound.

    Several upstream identifiers can land on one InChIKey. The label is defined as the
    majority call across them; a tie is dropped.
    """
    df = _standardise(records, structures)

    rows: list[dict] = []
    for inchikey, group in df.groupby("inchikey"):
        active_frac = float(group["active"].mean())
        if active_frac == 0.5:
            continue
        reported = group["potency_um"].dropna()
        rows.append(
            {
                "inchikey": inchikey,
                "smiles": group["smiles"].iloc[0],
                "label": int(active_frac > 0.5),
                # Provenance only: the label is the assay's own call, not a
                # threshold on this column.
                "potency_um": float(np.median(reported)) if len(reported) else float("nan"),
                "n_calls": len(group),
                "active_frac": round(active_frac, 6),
            }
        )
    return pd.DataFrame(rows).sort_values("inchikey").reset_index(drop=True)


def _counter_screen(records: pd.DataFrame, structures: pd.DataFrame) -> pd.DataFrame:
    """Collapse viability records the same way, keyed for a join on InChIKey."""
    df = _standardise(records, structures)

    rows: list[dict] = []
    for inchikey, group in df.groupby("inchikey"):
        active_frac = float(group["active"].mean())
        if active_frac == 0.5:
            continue
        rows.append(
            {
                "inchikey": inchikey,
                "cytotox": int(active_frac > 0.5),
                "cytotox_n_calls": len(group),
                "cytotox_active_frac": round(active_frac, 6),
            }
        )
    return pd.DataFrame(rows)


def _pool_table(records: dict, structures: pd.DataFrame) -> pd.DataFrame:
    """One row per compound, one call column per further viability screen."""
    from vp_core import dataset

    merged: pd.DataFrame | None = None
    smiles: dict[str, str] = {}
    for endpoint, frame in records.items():
        standardised = _standardise(frame, structures)
        smiles.update(
            dict(zip(standardised["inchikey"], standardised["smiles"], strict=True))
        )
        grouped = standardised.groupby("inchikey")["active"].mean()
        # The same rule the primary table uses: majority call, a tie dropped.
        called = grouped[grouped != 0.5]
        column = pd.DataFrame(
            {
                "inchikey": called.index,
                f"y_{endpoint.name.lower()}": (called > 0.5).astype(int).to_numpy(),
            }
        )
        merged = column if merged is None else merged.merge(column, on="inchikey", how="outer")

    if merged is None:
        raise RuntimeError("no viability screens returned any record")
    merged["smiles"] = merged["inchikey"].map(smiles)
    merged = merged.dropna(subset=["smiles"])
    for name in POOL_LABELS:
        merged[name] = merged[name].astype("Int64")
    ordered = ["inchikey", "smiles", *POOL_LABELS]
    table = merged[ordered].sort_values("inchikey").reset_index(drop=True)

    problems = dataset.validate_table(table, labels=POOL_LABELS)
    if problems:
        raise ValueError(f"rebuilt viability pool is invalid: {'; '.join(problems)}")
    return table


def fetch(*, write: bool = True) -> pd.DataFrame:
    """Rebuild the dataset from PubChem. Returns the standardised table."""
    from vp_core import dataset

    print(f"fetching AID {TARGET.pubchem_aid} from PubChem", file=sys.stderr)
    primary = _assay_records(TARGET.pubchem_aid)
    print(f"fetching AID {CYTOTOX.pubchem_aid} from PubChem", file=sys.stderr)
    counter = _assay_records(CYTOTOX.pubchem_aid)

    pool_records = {}
    for endpoint in VIABILITY_POOL:
        print(f"fetching AID {endpoint.pubchem_aid} from PubChem", file=sys.stderr)
        pool_records[endpoint] = _assay_records(endpoint.pubchem_aid)

    # One structure resolution for every assay, so the SMILES behind a
    # compound is the same string wherever it appears.
    cids = sorted(
        set(primary["cid"])
        | set(counter["cid"])
        | {c for frame in pool_records.values() for c in frame["cid"]}
    )
    structures = _smiles_for(cids)

    table = _to_compounds(primary, structures)
    # A left join: the reporter endpoint decides which compounds the table
    # holds, and the counter-screen fills in where it also has a call.
    table = table.merge(_counter_screen(counter, structures), on="inchikey", how="left")
    table["cytotox"] = table["cytotox"].astype("Int64")

    problems = dataset.validate_table(table, labels=LABELS)
    if problems:
        raise ValueError(f"rebuilt table is invalid: {'; '.join(problems)}")
    if write:
        dataset.write_table(table, TABLE_PATH, labels=LABELS)

    pool = _pool_table(pool_records, structures)
    if write:
        dataset.write_table(pool, POOL_PATH, labels=POOL_LABELS)
    print(
        f"viability pool: {len(pool)} compounds over {len(POOL_LABELS)} screens\n"
        f"  sha256 {dataset.dataset_hash(pool, labels=POOL_LABELS)}",
        file=sys.stderr,
    )

    called = table["cytotox"].notna()
    print(
        f"{len(table)} compounds, {int(table['label'].sum())} positive "
        f"({table['label'].mean():.1%})\n"
        f"  counter-screen calls {int(called.sum())} of them, "
        f"{int(table.loc[called, 'cytotox'].sum())} positive\n"
        f"  sha256 {dataset.dataset_hash(table, labels=LABELS)}",
        file=sys.stderr,
    )
    return table


def verify(version: str | None = None) -> dict:
    """Compare the on-disk table against the hash a released version recorded."""
    import vp_nrf2
    from vp_core import dataset
    from vp_core import manifest as manifest_mod

    resolved = vp_nrf2.get(version)
    labels = manifest_mod.dataset_labels(resolved.manifest)
    declared = resolved.manifest.get("dataset", {}).get("sha256")
    actual = dataset.dataset_hash(load(), labels=labels)
    out = {
        "version": resolved.name,
        "labels": labels,
        "declared": declared,
        "actual": actual,
        "match": declared == actual,
    }

    auxiliary = []
    for entry in resolved.manifest.get("dataset", {}).get("auxiliary", []):
        entry_labels = list(entry.get("labels", POOL_LABELS))
        path = entry.get("path")
        if path == "data/nrf2_epa_aeid1110.parquet":
            table = load_epa()
        elif path == "data/nrf2_viability_pool.parquet":
            table = load_pool()
        else:
            raise ValueError(f"unknown NRF2 auxiliary dataset: {path}")
        found = dataset.dataset_hash(table, labels=entry_labels)
        auxiliary.append(
            {
                "name": entry.get("name"),
                "labels": entry_labels,
                "declared": entry.get("sha256"),
                "actual": found,
                "match": entry.get("sha256") == found,
            }
        )
    if auxiliary:
        out["auxiliary"] = auxiliary
        out["match"] = out["match"] and all(a["match"] for a in auxiliary)
    return out


def build_example(n: int = 200, seed: int = 0) -> pd.DataFrame:
    """Regenerate the committed fixture from the full table."""
    from vp_core import dataset

    sample = dataset.stratified_example(load(), n=n, seed=seed, labels=LABELS)
    dataset.write_table(sample, EXAMPLE_PATH, labels=LABELS)

    # The fixture for the pool holds the same compounds the primary fixture
    # holds, so a test fitting the pooled head runs offline on both.
    pool = load_pool()
    held = list(sample["inchikey"])
    keep = pool[pool["inchikey"].isin(held)]
    extra = pool[~pool["inchikey"].isin(held)].head(len(sample))
    subset = (
        pd.concat([keep, extra], ignore_index=True)
        .sort_values("inchikey")
        .reset_index(drop=True)
    )
    dataset.write_table(subset, POOL_EXAMPLE_PATH, labels=POOL_LABELS)
    epa = dataset.stratified_example(load_epa(), n=n, seed=seed, labels=("label",))
    dataset.write_table(epa, EPA_EXAMPLE_PATH, labels=("label",))
    return sample


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m vp_nrf2.data")
    sub = parser.add_subparsers(dest="command", required=True)

    fetch_cmd = sub.add_parser("fetch", help="rebuild the Tox21 reporter and viability tables")
    fetch_cmd.add_argument(
        "--verify",
        action="store_true",
        help="after fetching, compare the hash against the released version",
    )
    sub.add_parser("verify", help="check the on-disk table against the recorded hash")
    epa_cmd = sub.add_parser("fetch-epa", help="rebuild the EPA ARE training table")
    epa_cmd.add_argument("--verify", action="store_true")
    example_cmd = sub.add_parser("build-example", help="regenerate the test fixture")
    example_cmd.add_argument("--n", type=int, default=200)

    args = parser.parse_args(argv)

    if args.command == "fetch":
        fetch()
        if not args.verify:
            return 0
    if args.command == "fetch-epa":
        from vp_nrf2 import epa

        epa.fetch()
        if not args.verify:
            return 0
    if args.command == "build-example":
        sample = build_example(n=args.n)
        print(f"wrote {len(sample)} rows to {EXAMPLE_PATH}")
        return 0

    report = verify()
    if report["match"]:
        print(f"dataset matches {report['version']}: {report['actual']}")
        return 0
    print(
        f"DRIFT: {report['version']} recorded {report['declared']}\n"
        f"       the rebuilt table hashes to {report['actual']}\n"
        "The upstream source has changed. Record it as a new version rather than "
        "overwriting the released one.",
        file=sys.stderr,
    )
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
