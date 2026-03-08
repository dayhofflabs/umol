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
| H | yes | yes (+1) | yes | yes (`atoms/h.toml`, `ions/h+1.toml`, hydrides) | review |  |
| He | yes | no | n/a | yes (`atoms/he.toml`) | review | |
| Li | yes | yes (+1) | no | yes (`atoms/li.toml`, `ions/li+1.toml`) | review |  |
| Be | yes | yes (+2) | no | yes (`atoms/be.toml`, `ions/be+2.toml`) | review |  |
| B | yes | no | no | yes (`atoms/b.toml`) | review |  |
| C | yes | no | no | yes (`atoms/c.toml`) | review |  |
| N | yes | yes (-3) | no | yes (`atoms/n.toml`, `ions/n-3.toml`) | review |  |
| O | yes | yes (-2) | no | yes (`atoms/o.toml`, `ions/o-2.toml`) | review |  |
| F | yes | yes (-1) | no | yes (`atoms/f.toml`, `ions/f-1.toml`, hydrides) | review |  |
| Ne | yes | no | n/a | yes (`atoms/ne.toml`) | review |  |
| Na | yes | yes (+1) | no | yes (`atoms/na.toml`, `ions/na+1.toml`) | review |  |
| Mg | yes | yes (+2) | no | yes (`atoms/mg.toml`, `ions/mg+2.toml`) | review |  |
| Al | yes | yes (+3) | no | yes (`atoms/al.toml`, `ions/al+3.toml`) | review |  |
| Si | yes | no | no | yes (`atoms/si.toml`) | review |  |
| P | yes | yes (-3) | no | yes (`atoms/p.toml`, `ions/p-3.toml`) | review |  |
| S | yes | yes (-2) | no | yes (`atoms/s.toml`, `ions/s-2.toml`) | review |  |
| Cl | yes | yes (-1) | no | yes (`atoms/cl.toml`, `ions/cl-1.toml`) | review |  |
| Ar | yes | no | n/a | yes (`atoms/ar.toml`) | review |  |
| K | yes | yes (+1) | no | yes (`atoms/k.toml`, `ions/k+1.toml`) | review |  |
| Ca | yes | yes (+2) | no | yes (`atoms/ca.toml`, `ions/ca+2.toml`) | review |  |
| Sc | yes | yes (+3) | no | yes (`atoms/sc.toml`, `ions/sc+3.toml`) | review |  |
| Ti | yes | yes (+3,+4) | no | yes (`atoms/ti.toml`, `ions/ti+3.toml`, `ions/ti+4.toml`) | review |  |
| V | yes | yes (+2,+3) | no | yes (`atoms/v.toml`, `ions/v+2.toml`, `ions/v+3.toml`) | review | Higher oxidation states intentionally excluded from ion set |
| Cr | yes | yes (+2,+3) | no | yes (`atoms/cr.toml`, `ions/cr+2.toml`, `ions/cr+3.toml`) | review | Higher oxidation states intentionally excluded from ion set |
| Mn | yes | yes (+2,+3) | no | yes (`atoms/mn.toml`, `ions/mn+2.toml`, `ions/mn+3.toml`) | review | Higher oxidation states intentionally excluded from ion set |
| Fe | yes | yes (+2,+3) | no | yes (`atoms/fe.toml`, `ions/fe+2.toml`, `ions/fe+3.toml`) | review |  |
| Co | yes | yes (+2,+3) | no | yes (`atoms/co.toml`, `ions/co+2.toml`, `ions/co+3.toml`) | review |  |
| Ni | yes | yes (+2) | no | yes (`atoms/ni.toml`, `ions/ni+2.toml`) | review |  |
| Cu | yes | yes (+1,+2) | no | yes (`atoms/cu.toml`, `ions/cu+1.toml`, `ions/cu+2.toml`) | review |  |
| Zn | yes | yes (+2) | no | yes (`atoms/zn.toml`, `ions/zn+2.toml`) | review |  |
| Ga | yes | yes (+3) | no | yes (`atoms/ga.toml`, `ions/ga+3.toml`) | review |  |
| Ge | yes | no | no | yes (`atoms/ge.toml`) | review |  |
| As | yes | yes (-3) | no | yes (`atoms/as.toml`, `ions/as-3.toml`) | review |  |
| Se | yes | yes (-2) | no | yes (`atoms/se.toml`, `ions/se-2.toml`) | review |  |
| Br | yes | yes (-1) | no | yes (`atoms/br.toml`, `ions/br-1.toml`, hydrides) | review |  |
| Kr | yes | no | n/a | yes (`atoms/kr.toml`) | review |  |
| Rb | yes | yes (+1) | no | yes (`atoms/rb.toml`, `ions/rb+1.toml`) | review |  |
| Sr | yes | yes (+2) | no | yes (`atoms/sr.toml`, `ions/sr+2.toml`) | review |  |
| Y | yes | yes (+3) | no | yes (`atoms/y.toml`, `ions/y+3.toml`) | review |  |
| Zr | yes | yes (+4) | no | yes (`atoms/zr.toml`, `ions/zr+4.toml`) | review |  |
| Nb | yes | yes (+3) | no | yes (`atoms/nb.toml`, `ions/nb+3.toml`) | review |  |
| Mo | yes | yes (+2,+3) | no | yes (`atoms/mo.toml`, `ions/mo+2.toml`, `ions/mo+3.toml`) | review |  |
| Tc | yes | no | no | yes (`atoms/tc.toml`) | review | No ion added in current scope |
| Ru | yes | yes (+2,+3) | no | yes (`atoms/ru.toml`, `ions/ru+2.toml`, `ions/ru+3.toml`) | review |  |
| Rh | yes | yes (+2,+3) | no | yes (`atoms/rh.toml`, `ions/rh+2.toml`, `ions/rh+3.toml`) | review |  |
| Pd | yes | yes (+2) | no | yes (`atoms/pd.toml`, `ions/pd+2.toml`) | review |  |
| Ag | yes | yes (+1,+2) | no | yes (`atoms/ag.toml`, `ions/ag+1.toml`, `ions/ag+2.toml`) | review |  |
| Cd | yes | yes (+2) | no | yes (`atoms/cd.toml`, `ions/cd+2.toml`) | review |  |
| In | yes | yes (+3) | no | yes (`atoms/in.toml`, `ions/in+3.toml`) | review |  |
| Sn | yes | no | no | yes (`atoms/sn.toml`) | review |  |
| Sb | yes | yes (-3) | no | yes (`atoms/sb.toml`, `ions/sb-3.toml`) | review |  |
| Te | yes | yes (-2) | no | yes (`atoms/te.toml`, `ions/te-2.toml`) | review |  |
| I | yes | yes (-1) | no | yes (`atoms/i.toml`, `ions/i-1.toml`, hydrides) | review | `implicit_h=false => H=Some(0)` rule required for unambiguous typing |
| Xe | yes | no | partial | yes (`atoms/xe.toml`) | review |  |

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

## Current Conformance Set

- `atoms/`: `h`, `he`, `li`, `be`, `b`, `c`, `n`, `o`, `f`, `ne`, `na`, `mg`, `al`, `si`, `p`, `s`, `cl`, `ar`, `k`, `ca`, `sc`, `ti`, `v`, `cr`, `mn`, `fe`, `co`, `ni`, `cu`, `zn`, `ga`, `ge`, `as`, `se`, `br`, `kr`, `rb`, `sr`, `y`, `zr`, `nb`, `mo`, `tc`, `ru`, `rh`, `pd`, `ag`, `cd`, `in`, `sn`, `sb`, `te`, `i`, `xe`
- `ions/`: `h+1`, `li+1`, `be+2`, `mg+2`, `al+3`, `n-3`, `o-2`, `f-1`, `na+1`, `p-3`, `s-2`, `cl-1`, `k+1`, `ca+2`, `sc+3`, `ti+3`, `ti+4`, `v+2`, `v+3`, `cr+2`, `cr+3`, `mn+2`, `mn+3`, `fe+2`, `fe+3`, `co+2`, `co+3`, `ni+2`, `cu+1`, `cu+2`, `zn+2`, `ga+3`, `as-3`, `se-2`, `br-1`, `rb+1`, `sr+2`, `y+3`, `zr+4`, `nb+3`, `mo+2`, `mo+3`, `ru+2`, `ru+3`, `rh+2`, `rh+3`, `pd+2`, `ag+1`, `ag+2`, `cd+2`, `in+3`, `sb-3`, `te-2`, `i-1`
- `hydrides/`: `h2`, `hf`

## Decision Log

### 2026-03-07

- Agreed staged rollout: ground states -> discrete ions -> typical non-metal valences.
- Agreed registry expansion must be paired with conformance additions.
- Agreed registry specs are reviewed with code-level rigor (small, explicit, traceable changes).
- Agreed conformance naming/categories: `atoms/`, `ions/`, `hydrides/` and file naming rules above.
- Added `implicit_h` control to conformance inputs; counts strategy now honors this via `enable_implicit_hydrogens`.
- Expanded coverage to period-3 atoms (`Na`..`Ar`) and matching ions (`Na+`, `Al3+`, `P3-`, `S2-`, `Cl-`), plus `Be2+`.
- Expanded coverage to period-4 and period-5 s/p blocks, 3d block, and 4d block with ions: `Y3+`, `Zr4+`, `Nb3+`, `Mo2+`, `Mo3+`, `Ru2+`, `Ru3+`, `Rh2+`, `Rh3+`, `Pd2+`, `Ag+`, `Ag2+`, `Cd2+`.
- Added missing 3d atom conformance files (`atoms/sc.toml` through `atoms/zn.toml`) and refreshed snapshots.
- TODO: 5d ions La3+, Hf4+, Au3+, Hg2+, [Hg2]2+ ("{Hg+/5v1}"), 