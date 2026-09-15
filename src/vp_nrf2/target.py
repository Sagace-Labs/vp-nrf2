"""The endpoints this pathway predicts.

Frozen so the data sources are auditable.

Both assays were verified live against PubChem on 2026-09-06, and screen the
same library:

    AID 743202   qHTS assay for small molecule agonists of the antioxidant
                 response element (ARE) signaling pathway  9305 substances
    AID 743203   the same, cell viability counter screen    9305 substances
"""

from __future__ import annotations

from dataclasses import dataclass

__all__ = ["CYTOTOX", "TARGET", "TARGETS", "Endpoint", "all_names", "get"]


@dataclass(frozen=True)
class Endpoint:
    """One molecular initiating event, identified by the assay that reads it out."""

    name: str
    pathway: str
    mie: str
    pubchem_aid: int
    assay_name: str


TARGET = Endpoint(
    name="ARE",
    pathway="NRF2 / KEAP1 antioxidant response",
    mie=(
        "covalent modification of reactive KEAP1 cysteines by an electrophile, "
        "which releases NRF2 and switches on the antioxidant response element"
    ),
    pubchem_aid=743202,
    assay_name=(
        "qHTS assay for small molecule agonists of the antioxidant response "
        "element (ARE) signaling pathway"
    ),
)

CYTOTOX = Endpoint(
    name="VIABILITY",
    pathway="NRF2 / KEAP1 antioxidant response",
    mie="loss of cell viability, which registers on the reporter readout as a consequence",
    pubchem_aid=743203,
    assay_name=(
        "qHTS assay for small molecule agonists of the antioxidant response "
        "element (ARE) signaling pathway - cell viability counter screen"
    ),
)

TARGETS: dict[str, Endpoint] = {"ARE": TARGET, "VIABILITY": CYTOTOX}


def get(name: str) -> Endpoint:
    try:
        return TARGETS[name.upper()]
    except KeyError:
        raise KeyError(f"unknown target {name!r}; known: {sorted(TARGETS)}") from None


def all_names() -> list[str]:
    return sorted(TARGETS)
