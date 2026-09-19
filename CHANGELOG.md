# Changelog

Released versions are immutable. A correction to a released version is a new
patch version. This file records *why* each version exists; what
it contains is in its own `manifest.toml` and `CARD.md`.

## v2

`nrf2_cytotox` trains on the panel's five other Tox21 viability
counter-screens as well as its own. Its label is unchanged: the endpoint is
still the call AID 743203 makes. Its own screen supplies 40 of the positives
the head had to learn from, which is too few.

Donor rows are blocked by Murcko scaffold against the fold that judges the
fit, so a held-out compound and its analogues are excluded even where another
screen labels them.

Signature 1, unchanged.

## v1

Initial release.

Signature 1: `nrf2_are` and `nrf2_cytotox`.
