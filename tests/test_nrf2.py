"""NRF2 package tests.

Everything here runs on the committed fixture and the shipped weights.
"""

from __future__ import annotations

import numpy as np
import pytest

import vp_nrf2
from vp_nrf2 import contract, data

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
CAFFEINE = "Cn1cnc2c1c(=O)n(C)c(=O)n2C"
NONSENSE = "not-a-molecule"


def test_versions_are_discoverable_and_ordered():
    versions = vp_nrf2.versions()
    assert versions, "no released versions found"
    assert vp_nrf2.current_version() == versions[-1]


def test_declared_signature_matches_the_shipped_manifest():
    shipped = vp_nrf2.signature()
    assert shipped["version"] == contract.SIGNATURE_VERSION
    assert [o["name"] for o in shipped["outputs"]] == contract.column_names()


def test_predictions_carry_the_declared_column_names():
    frame = vp_nrf2.predict([ASPIRIN, CAFFEINE])
    assert list(frame.columns) == contract.column_names()
    assert len(frame) == 2


def test_probabilities_are_in_range():
    frame = vp_nrf2.predict([ASPIRIN, CAFFEINE])
    for column in contract.column_names():
        values = frame[column].to_numpy()
        assert np.all((values >= 0.0) & (values <= 1.0))


def test_unparseable_input_becomes_nan_rather_than_raising():
    frame = vp_nrf2.predict([ASPIRIN, NONSENSE])
    assert np.isfinite(frame[contract.PRIMARY].iloc[0])
    assert frame.iloc[1].isna().all()


def test_a_version_can_be_pinned():
    """Every release stays loadable under its own contract, not the newest one."""
    for name in vp_nrf2.versions():
        frame = vp_nrf2.predict([ASPIRIN], version=name)
        declared = [o["name"] for o in vp_nrf2.signature(name)["outputs"]]
        assert list(frame.columns) == declared

    assert vp_nrf2.predict([ASPIRIN]).equals(
        vp_nrf2.predict([ASPIRIN], version=vp_nrf2.current_version())
    )


def test_the_two_endpoints_are_distinct_predictions():
    """A second column that merely copied the first would buy nothing."""
    frame = vp_nrf2.predict([ASPIRIN, CAFFEINE, "CCCCCCCCCCCCn1cc[n+](C)c1"])
    assert not np.allclose(
        frame["nrf2_are"].to_numpy(), frame["nrf2_cytotox"].to_numpy()
    )


def test_a_single_string_is_rejected():
    with pytest.raises(TypeError, match="list of SMILES"):
        vp_nrf2.predict(ASPIRIN)  # type: ignore[arg-type]


def test_unknown_version_names_the_available_ones():
    with pytest.raises(FileNotFoundError, match="v1"):
        vp_nrf2.predict([ASPIRIN], version="v999")


def test_fixture_satisfies_the_dataset_contract():
    from vp_core import dataset as core_dataset

    fixture = data.example()
    assert core_dataset.validate_table(fixture, labels=data.LABELS) == []
    for column in data.LABELS:
        assert set(fixture[column].dropna()) == {0, 1}


def test_an_uncalled_endpoint_is_absent_rather_than_negative():
    """The distinction the second label column exists to preserve."""
    fixture = data.example()
    assert fixture["cytotox"].isna().any()
    called = data.labelled(fixture, "cytotox")
    assert len(called) == int(fixture["cytotox"].notna().sum())
    assert len(called) < len(fixture)


def test_majority_vote_drops_a_compound_the_assay_called_both_ways():
    """A tie carries no clean label, so it must not be broken arbitrarily."""
    import pandas as pd

    records = pd.DataFrame(
        {
            "cid": [1, 2, 3, 3],
            "outcome": ["Active", "Inactive", "Active", "Inactive"],
            "potency_um": [4.0, float("nan"), 8.0, float("nan")],
        }
    )
    structures = pd.DataFrame(
        {"cid": [1, 2, 3], "smiles_raw": [ASPIRIN, CAFFEINE, "CCO"]}
    )
    table = data._to_compounds(records, structures)

    assert len(table) == 2
    assert set(table["label"]) == {0, 1}
    assert "CCO" not in set(table["smiles"])


def test_target_registry_is_uniform():
    assert vp_nrf2.all_names() == ["ARE", "VIABILITY"]
    assert vp_nrf2.get_target("are").pubchem_aid == 743202
    assert vp_nrf2.get_target("viability").pubchem_aid == 743203
    with pytest.raises(KeyError):
        vp_nrf2.get_target("nope")


def test_evaluation_harness_runs_on_the_fixture():
    """A smoke run must never be recordable as a released result."""
    from vp_nrf2.evaluate import evaluate_version

    try:
        record = evaluate_version(vp_nrf2.current_version(), use_example=True, write=False)
    except ValueError as exc:
        if "produced an empty" not in str(exc):
            raise
        pytest.skip(f"the fixture cannot be split under this protocol: {exc}")
    assert record["dataset"] == "example fixture"
    assert [o["name"] for o in record["additional_outputs"]] == [
        name for name in contract.column_names() if name != contract.PRIMARY
    ]
    with pytest.raises(ValueError, match="refusing to record"):
        evaluate_version(vp_nrf2.current_version(), use_example=True, write=True)
