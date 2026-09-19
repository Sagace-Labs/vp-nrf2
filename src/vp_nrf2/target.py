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

__all__ = [
    "CYTOTOX",
    "TARGET",
    "TARGETS",
    "VIABILITY_POOL",
    "Endpoint",
    "all_names",
    "get",
]


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

#: Further Tox21 viability counter-screens, fetched by this package and used
#: as extra training rows for ``nrf2_cytotox``. They are listed here as assay
#: ids rather than imported from the packages that also read them.
#: Verified live 2026-09-06.
VIABILITY_POOL: tuple[Endpoint, ...] = (
    Endpoint(
        name="AHR_VIABILITY",
        pathway="aryl hydrocarbon receptor",
        mie="loss of cell viability",
        pubchem_aid=743086,
        assay_name=(
            "qHTS assay to identify small molecule that activate the aryl "
            "hydrocarbon receptor (AhR) signaling pathway - cell viability "
            "counter screen"
        ),
    ),
    Endpoint(
        name="GR_VIABILITY",
        pathway="glucocorticoid receptor",
        mie="loss of cell viability",
        pubchem_aid=720693,
        assay_name=(
            "qHTS assay to identify small molecule antagonists of the "
            "glucocorticoid receptor (GR) signaling pathway - cell viability "
            "counter screen"
        ),
    ),
    Endpoint(
        name="HSR_VIABILITY",
        pathway="heat shock response",
        mie="loss of cell viability",
        pubchem_aid=743209,
        assay_name=(
            "qHTS assay for small molecule activators of the heat shock "
            "response signaling pathway - cell viability counter screen"
        ),
    ),
    Endpoint(
        name="P53_VIABILITY",
        pathway="p53 genotoxic stress response",
        mie="loss of cell viability",
        pubchem_aid=651633,
        assay_name=(
            "qHTS assay to identify small molecule agonists of the p53 "
            "signaling pathway - cell viability counter screen"
        ),
    ),
    Endpoint(
        name="PXR_VIABILITY",
        pathway="pregnane X receptor",
        mie="loss of cell viability",
        pubchem_aid=1346977,
        assay_name=(
            "qHTS assay to identify small molecule agonists of the pregnane X "
            "receptor (PXR) signaling pathway - cell viability counter screen"
        ),
    ),
)


def get(name: str) -> Endpoint:
    try:
        return TARGETS[name.upper()]
    except KeyError:
        raise KeyError(f"unknown target {name!r}; known: {sorted(TARGETS)}") from None


def all_names() -> list[str]:
    return sorted(TARGETS)
