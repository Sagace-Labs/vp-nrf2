# NRF2 data

    data/
      nrf2_tox21.parquet                    standardised table (producing the hash)
      nrf2_viability_pool.parquet           further viability screens, training only
      example/
        nrf2_example.parquet                small stratified fixture, used by the tests
        nrf2_viability_pool_example.parquet the same for the pool

`nrf2_tox21.parquet` holds one row per compound with the two identity columns —
`inchikey` and `smiles` (standardised) — and two label columns. `label` is the
antioxidant-response reporter call and `cytotox` is the viability call.
Alongside them are associated `potency_um`, `n_calls`, `active_frac`,
`cytotox_n_calls` and `cytotox_active_frac` as provenance. Only the identity
and label columns are hashed, so a provenance column may be added without
moving the dataset hash.

## Origin and processing

PubChem BioAssay AID 743202, the Tox21 qHTS screen for agonists of the
antioxidant response element signalling pathway, and AID 743203, the
cell-viability counter screen over the same library, retrieved through the
public PUG-REST concise assay endpoint. Rows each assay called Active or
Inactive are kept and rows it called Inconclusive are dropped. Compound
identifiers are resolved to isomeric SMILES through the compound property
endpoint; structures are standardised (normalise, largest fragment,
neutralise); records are collapsed to one row per InChIKey by majority call,
and a tie is dropped.

A standardised structure that RDKit will not read back is dropped too, in order
to avoid featurization outputting all-zeros.

The reporter screen decides which compounds the table holds. A compound the
counter-screen did not call carries a null in `cytotox`, which records that no
call was made.

`nrf2_viability_pool.parquet` holds the viability counter-screens of five
further Tox21 screens — AIDs 743086, 720693, 743209, 651633 and 1346977 —
processed the same way, one call column per screen and a null where that
screen made no call. It supplies extra training rows for `nrf2_cytotox` and is
never evaluated against; the endpoint remains AID 743203's call. It carries
its own hash in the manifest under `[[dataset.auxiliary]]`.

`potency_um` is the median of the potencies reported across a compound's
reporter records and is NaN when none reported one. `n_calls` is how many assay
records collapsed into the row and `active_frac` is the share of them that were
Active, so the majority vote stays auditable. The two `cytotox_` columns say the
same of the counter-screen.

## Rebuilding it

    python -m vp_nrf2.data fetch --verify

This downloads, re-parses and re-hashes, then compares against the hash the
current version recorded. A mismatch in the hash can indicate upstream changes
to the source data or the processing pipeline.

## Licence

The bundled table is a United States government work in the public domain; see
`../LICENSE-DATA`. That file covers these files only, not the package code or
the trained weights.
