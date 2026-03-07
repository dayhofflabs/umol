# Default Registry Rollout Plan (2026-03-07)

This document tracks incremental, reviewable expansion of the default atom type registry.
Registry entries are treated as code: explicit scope, per-entry review, and matching conformance coverage.

## Scope and Ordering

Implementation proceeds in this order:

1. Atomic ground states
2. Well-defined discrete ions (especially metals)
3. Typical valences for non-metals

## Review Rules

- Every new registry spec must be accompanied by conformance input(s).
- No bulk "periodic table dump" commits.
- Each element is added in a small, auditable change set.
- Deviations from expected chemistry behavior are documented explicitly.
- Registry edits and conformance edits are reviewed together.

## Change Unit (Per Element)

For each element `E`, a minimal unit of work is:

1. Add/adjust `E` entries in `umol-models-graph/config/default-registry.toml`
2. Add/adjust `E` conformance files in `umol-models-graph/tests/resolution/data/`
3. Update/accept corresponding snapshots
4. Record status in this document

## Tracking Table

Status values:

- `todo`: not started
- `in_progress`: currently being edited
- `review`: implemented, awaiting chemistry review
- `done`: accepted

| Element | Ground state | Discrete ions | Typical valences | Conformance added | Status | Notes |
|---|---|---|---|---|---|---|
| H | yes | yes (+1) | yes | yes (`atoms/h.toml`, `ions/h+1.toml`, hydrides) | review | Counts now respects `implicit_h`; verify final ground-state convention |
| He | yes | no | n/a | yes (`atoms/he.toml`) | review | Counts currently yields `{He/1}` |
| Li | yes | yes (+1) | no | yes (`atoms/li.toml`, `ions/li+1.toml`) | review | Counts/typing parity fixed under `implicit_h=false` |
| Be | yes | yes (+2) | no | yes (`atoms/be.toml`, `ions/be+2.toml`) | review | Added missing `Be+2` registry spec |
| B | yes | no | no | yes (`atoms/b.toml`) | review |  |
| C | yes | no | no | yes (`atoms/c.toml`) | review |  |
| N | yes | yes (-3) | no | yes (`atoms/n.toml`, `ions/n-3.toml`) | review | `N3-` included intentionally (hypothetical edge case) |
| O | yes | yes (-2) | no | yes (`atoms/o.toml`, `ions/o-2.toml`) | review | `O2-` included intentionally (hypothetical edge case) |
| F | yes | yes (-1) | no | yes (`atoms/f.toml`, `ions/f-1.toml`, hydrides) | review |  |
| Ne | yes | no | n/a | yes (`atoms/ne.toml`) | review |  |
| Na | yes | yes (+1) | no | yes (`atoms/na.toml`, `ions/na+1.toml`) | review |  |
| Mg | yes | yes (+2) | no | yes (`atoms/mg.toml`, `ions/mg+2.toml`) | review |  |
| Al | yes | yes (+3) | no | yes (`atoms/al.toml`, `ions/al+3.toml`) | review |  |
| Si | yes | no | no | yes (`atoms/si.toml`) | review |  |
| P | yes | partial (-3 queried) | no | yes (`atoms/p.toml`, `ions/p-3.toml`) | review | `P3-` missing in atom-typing registry (counts-only pass) |
| S | yes | yes (-2) | no | yes (`atoms/s.toml`, `ions/s-2.toml`) | review |  |
| Cl | yes | yes (-1) | no | yes (`atoms/cl.toml`, `ions/cl-1.toml`, hydrides) | review |  |
| Ar | yes | no | n/a | yes (`atoms/ar.toml`) | review |  |

Extend this table as new elements are considered.

## Conformance Policy

Per element, add explicit resolution inputs that cover:

- Free atom query (with explicit `implicit_h` policy)
- At least one bonded case when chemically meaningful
- Relevant charge states (for ion stage)
- At least one "should fail" case where useful to constrain over-broad matching

Suggested naming pattern:

- Atomic ground state: `tests/resolution/data/atoms/<element>.toml`
- Atomic valence/excited variants: `tests/resolution/data/atoms/<element>_<motif>.toml`
- Ionic ground state: `tests/resolution/data/ions/<element><charge>.toml`
- Ionic valence/excited variants: `tests/resolution/data/ions/<element><charge>_<motif>.toml`

Conventions:

- `<element>`: lowercase element symbol/name token used in file naming.
- `<charge>`: signed integer with explicit sign, e.g. `+1`, `-2` (no underscore between element and charge).
- `<motif>`: lowercase, underscore-separated descriptor; may encode occupation and optional multiplet.

## Categories

Use resolution data categories:

- `atoms/`
- `ions/`
- `hydrides/` (for existing molecular hydride cases moved from `basic/`)

## Current Conformance Seed

- `atoms/`: `h`, `he`, `li`, `be`, `b`, `c`, `n`, `o`, `f`, `ne`, `na`, `mg`, `al`, `si`, `p`, `s`, `cl`, `ar`
- `ions/`: `h+1`, `li+1`, `be+2`, `mg+2`, `al+3`, `n-3`, `o-2`, `f-1`, `na+1`, `p-3`, `s-2`, `cl-1`
- `hydrides/`: `h2`, `hf`

## Decision Log

### 2026-03-07

- Agreed staged rollout: ground states -> discrete ions -> typical non-metal valences.
- Agreed registry expansion must be paired with conformance additions.
- Agreed registry specs are reviewed with code-level rigor (small, explicit, traceable changes).
- Agreed conformance naming/categories: `atoms/`, `ions/`, `hydrides/` and file naming rules above.
- Added `implicit_h` control to conformance inputs; counts strategy now honors this via `enable_implicit_hydrogens`.
- Expanded coverage to period-3 atoms (`Na`..`Ar`) and matching ions (`Na+`, `Al3+`, `P3-`, `S2-`, `Cl-`), plus `Be2+`.
