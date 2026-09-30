# 229 — Python SMILES output

Status: Completed
Date: 2026-09-29
Relates: [153](153-format-parsing-outstanding-tasks-2026-07-18.md),
[170](170-reaction-smiles-python-2026-07-28.md),
[226](226-smiles-roundtrip-design-2026-09-10.md)

## Purpose

Doc 226 completed Rust SMILES output: `export_smiles`, `export_smiles_with`,
`export_reaction_smiles` and `export_reaction_smiles_with` in `umol-graph/src/export.rs`, with
stereo and, for reactions, the atom correspondence written as map labels. It left Python output to a
separate decision (doc 226 § Python, doc 153 T8). A downstream consumer, an editor that builds
reactions from `Reaction.from_sides` and hands them to a model that takes reaction SMILES, needs
that output from Python and has no other stereo-preserving writer.

## Design

Two methods, each the output counterpart of an existing ingestion method with the same keyword
arguments and the same defaults.

- `Molecule.to_smiles(self, *, io_config=None, chemistry_model=None, resolve_config=None) -> str`
  calls `export_smiles_with`.
- `Reaction.to_reaction_smiles(self, *, io_config=None, chemistry_model=None,
  resolve_config=None) -> str` snapshots the reaction and calls `export_reaction_smiles_with`.

Omitted options select what `Molecule.from_smiles` and `Reaction.from_reaction_smiles` select:
`SmilesIoConfig.opensmiles()`, the SMILES valence model, and the Natural isotope policy. These are
also the presets of Rust `export_smiles` / `export_reaction_smiles`, so default Python output is
the text the Rust convenience functions produce. No Python `Smiles` boundary value, `Convey`, or
`Resolver` is exposed; the methods return text only.

**Errors** follow the ingestion taxonomy. A projection contradiction, or a reaction that cannot be
materialized, raises `ContradictionError`; an underdetermined projection raises
`UnderdeterminedError`. Every other convey or render failure, where a determined value has no SMILES
representation (an unrepresentable field value, constraint, entity, stereo frame or map label, or a
render-unsupported field), raises `ModelConversionError`. The message is the Rust error's text,
including the `reactants:` / `products:` prefix for a side failure. The input is unchanged on
success and on failure.

## Consequences

The round-trip property of doc 226 is observable from Python: for the tested molecules and
reactions, `from_smiles(m.to_smiles())` and `from_reaction_smiles(r.to_reaction_smiles())` reproduce
the source text, including tetrahedral and cis/trans stereo. Map labels are renumbered from 1 in
correspondence order; original labels are not retained (doc 226).

Values built in Python must be determined to export. `Molecule.parse` leaves unwritten fields open,
so a parsed atom without an explicit isotope raises `ModelConversionError` ("non-ground isotope")
under the default policy. Molecules from `from_smiles` export as given.

## Tests

`umol-py/tests/test_smiles_export.py`: molecule round trips mirroring the Rust export cases; the
chemistry model reaching the exporter (the default `ChemistryModel` brackets every hydrogen count);
four `ModelConversionError` messages; reaction round trips including map-label renumbering and both
stereo kinds; a reaction built with `Reaction.from_sides`; a side-prefixed reaction error. Both
signatures are pinned in `tests/test_import.py`.
