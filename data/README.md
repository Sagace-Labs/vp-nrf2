# NRF2 data

    data/
      nrf2_tox21.parquet                    standardised table (producing the hash)
      nrf2_viability_pool.parquet           further viability screens, training only
      nrf2_epa_aeid1110.parquet             EPA v4.3 ARE calls, training only
      example/
        nrf2_example.parquet                small stratified fixture, used by the tests
        nrf2_viability_pool_example.parquet the same for the pool
        nrf2_epa_example.parquet            small stratified EPA training fixture

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

`nrf2_epa_aeid1110.parquet` contains 1,613 EPA invitrodb v4.3 AEID
1110 compounds absent from the primary table by connectivity identity. The
PubChem export of 18 August 2026 supplies Active/Inactive calls; Inconclusive
calls and unresolved structures are excluded. Structures come from EPA's
archived 2018 chemical identity SDF. Ties and conflicting stereoisomers are
dropped. During each fit, external rows sharing a scaffold with validation or
test compounds are also excluded. This table is training-only.

`potency_um` is the median of the potencies reported across a compound's
reporter records and is NaN when none reported one. `n_calls` is how many assay
records collapsed into the row and `active_frac` is the share of them that were
Active, so the majority vote stays auditable. The two `cytotox_` columns say the
same of the counter-screen.

## Rebuilding it

    python -m vp_nrf2.data fetch --verify
    python -m vp_nrf2.data fetch-epa --verify

This downloads, re-parses and re-hashes, then compares against the hash the
released version recorded. A mismatch in the hash can indicate upstream changes
to the source data or the processing pipeline.

## Licence

The bundled Tox21 and EPA tables carry no copyright restrictions; see
`../LICENSE-DATA`. That file covers data only, not the code or trained weights.
