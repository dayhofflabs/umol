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
| H | yes | yes (+1) | yes | yes (`atoms/h.toml`, `ions/h+1.toml`, hydrides) | done |  |
| He | yes | no | n/a | yes (`atoms/he.toml`) | done | |
| Li | yes | yes (+1) | no | yes (`atoms/li.toml`, `ions/li+1.toml`) | done |  |
| Be | yes | yes (+2) | no | yes (`atoms/be.toml`, `ions/be+2.toml`) | done |  |
| B | yes | no | no | yes (`atoms/b.toml`) | done |  |
| C | yes | no | no | yes (`atoms/c.toml`) | done |  |
| N | yes | yes (-3) | no | yes (`atoms/n.toml`, `ions/n-3.toml`) | done |  |
| O | yes | yes (-2) | no | yes (`atoms/o.toml`, `ions/o-2.toml`) | done |  |
| F | yes | yes (-1) | no | yes (`atoms/f.toml`, `ions/f-1.toml`, hydrides) | done |  |
| Ne | yes | no | n/a | yes (`atoms/ne.toml`) | done |  |
| Na | yes | yes (+1) | no | yes (`atoms/na.toml`, `ions/na+1.toml`) | done |  |
| Mg | yes | yes (+2) | no | yes (`atoms/mg.toml`, `ions/mg+2.toml`) | done |  |
| Al | yes | yes (+3) | no | yes (`atoms/al.toml`, `ions/al+3.toml`) | done |  |
| Si | yes | no | no | yes (`atoms/si.toml`) | done |  |
| P | yes | yes (-3) | no | yes (`atoms/p.toml`, `ions/p-3.toml`) | done |  |
| S | yes | yes (-2) | no | yes (`atoms/s.toml`, `ions/s-2.toml`) | done |  |
| Cl | yes | yes (-1) | no | yes (`atoms/cl.toml`, `ions/cl-1.toml`) | done |  |
| Ar | yes | no | n/a | yes (`atoms/ar.toml`) | done |  |
| K | yes | yes (+1) | no | yes (`atoms/k.toml`, `ions/k+1.toml`) | done |  |
| Ca | yes | yes (+2) | no | yes (`atoms/ca.toml`, `ions/ca+2.toml`) | done |  |
| Sc | yes | yes (+3) | no | yes (`atoms/sc.toml`, `ions/sc+3.toml`) | done |  |
| Ti | yes | yes (+3,+4) | no | yes (`atoms/ti.toml`, `ions/ti+3.toml`, `ions/ti+4.toml`) | done |  |
| V | yes | yes (+2,+3) | no | yes (`atoms/v.toml`, `ions/v+2.toml`, `ions/v+3.toml`) | done |  |
| Cr | yes | yes (+2,+3) | no | yes (`atoms/cr.toml`, `ions/cr+2.toml`, `ions/cr+3.toml`) | done |  |
| Mn | yes | yes (+2,+3) | no | yes (`atoms/mn.toml`, `ions/mn+2.toml`, `ions/mn+3.toml`) | done |  |
| Fe | yes | yes (+2,+3) | no | yes (`atoms/fe.toml`, `ions/fe+2.toml`, `ions/fe+3.toml`) | done |  |
| Co | yes | yes (+2,+3) | no | yes (`atoms/co.toml`, `ions/co+2.toml`, `ions/co+3.toml`) | done |  |
| Ni | yes | yes (+2) | no | yes (`atoms/ni.toml`, `ions/ni+2.toml`) | done |  |
| Cu | yes | yes (+1,+2) | no | yes (`atoms/cu.toml`, `ions/cu+1.toml`, `ions/cu+2.toml`) | done |  |
| Zn | yes | yes (+2) | no | yes (`atoms/zn.toml`, `ions/zn+2.toml`) | done |  |
| Ga | yes | yes (+3) | no | yes (`atoms/ga.toml`, `ions/ga+3.toml`) | done |  |
| Ge | yes | no | no | yes (`atoms/ge.toml`) | done |  |
| As | yes | yes (-3) | no | yes (`atoms/as.toml`, `ions/as-3.toml`) | done |  |
| Se | yes | yes (-2) | no | yes (`atoms/se.toml`, `ions/se-2.toml`) | done |  |
| Br | yes | yes (-1) | no | yes (`atoms/br.toml`, `ions/br-1.toml`, hydrides) | done |  |
| Kr | yes | no | n/a | yes (`atoms/kr.toml`) | done |  |
| Rb | yes | yes (+1) | no | yes (`atoms/rb.toml`, `ions/rb+1.toml`) | done |  |
| Sr | yes | yes (+2) | no | yes (`atoms/sr.toml`, `ions/sr+2.toml`) | done |  |
| Y | yes | yes (+3) | no | yes (`atoms/y.toml`, `ions/y+3.toml`) | done |  |
| Zr | yes | yes (+4) | no | yes (`atoms/zr.toml`, `ions/zr+4.toml`) | done |  |
| Nb | yes | yes (+3) | no | yes (`atoms/nb.toml`, `ions/nb+3.toml`) | done |  |
| Mo | yes | yes (+2,+3) | no | yes (`atoms/mo.toml`, `ions/mo+2.toml`, `ions/mo+3.toml`) | done |  |
| Tc | yes | no | no | yes (`atoms/tc.toml`) | done | No ion added in current scope |
| Ru | yes | yes (+2,+3) | no | yes (`atoms/ru.toml`, `ions/ru+2.toml`, `ions/ru+3.toml`) | done |  |
| Rh | yes | yes (+2,+3) | no | yes (`atoms/rh.toml`, `ions/rh+2.toml`, `ions/rh+3.toml`) | done |  |
| Pd | yes | yes (+2) | no | yes (`atoms/pd.toml`, `ions/pd+2.toml`) | done |  |
| Ag | yes | yes (+1,+2) | no | yes (`atoms/ag.toml`, `ions/ag+1.toml`, `ions/ag+2.toml`) | done |  |
| Cd | yes | yes (+2) | no | yes (`atoms/cd.toml`, `ions/cd+2.toml`) | done |  |
| In | yes | yes (+3) | no | yes (`atoms/in.toml`, `ions/in+3.toml`) | done |  |
| Sn | yes | no | no | yes (`atoms/sn.toml`) | done |  |
| Sb | yes | yes (-3) | no | yes (`atoms/sb.toml`, `ions/sb-3.toml`) | done |  |
| Te | yes | yes (-2) | no | yes (`atoms/te.toml`, `ions/te-2.toml`) | done |  |
| I | yes | yes (-1) | no | yes (`atoms/i.toml`, `ions/i-1.toml`, hydrides) | done |   |
| Xe | yes | no | partial | yes (`atoms/xe.toml`) | done |  |
| Cs | yes | yes (+1) | no | yes (`atoms/cs.toml`, `ions/cs+1.toml`) | done |  |
| Ba | yes | yes (+2) | no | yes (`atoms/ba.toml`, `ions/ba+2.toml`) | done |  |
| Tl | yes | yes (+1,+3) | no | yes (`atoms/tl.toml`, `ions/tl+1.toml`, `ions/tl+3.toml`) | done |  |
| Pb | yes | yes (+2,+4) | no | yes (`atoms/pb.toml`, `ions/pb+2.toml`, `ions/pb+4.toml`) | done |  |
| Bi | yes | yes (+3) | no | yes (`atoms/bi.toml`, `ions/bi+3.toml`) | done |  |
| Po | yes | yes (-2) | no | yes (`atoms/po.toml`, `ions/po-2.toml`) | done |  |
| At | yes | yes (-1) | no | yes (`atoms/at.toml`, `ions/at-1.toml`) | done |  |
| Rn | yes | no | no | yes (`atoms/rn.toml`) | done |  |
| La | yes | yes (+3) | no | yes (`atoms/la.toml`, `ions/la+3.toml`) | done |  |
| Ce | yes | yes (+3,+4) | no | yes (`atoms/ce.toml`, `ions/ce+3.toml`, `ions/ce+4.toml`) | done |  |
| Pr | yes | yes (+3,+4) | no | yes (`atoms/pr.toml`, `ions/pr+3.toml`, `ions/pr+4.toml`) | done | Counts requires `Pr` valence-table support for +4 |
| Nd | yes | yes (+3) | no | yes (`atoms/nd.toml`, `ions/nd+3.toml`) | done |  |
| Pm | yes | yes (+3) | no | yes (`atoms/pm.toml`, `ions/pm+3.toml`) | done |  |
| Sm | yes | yes (+2,+3) | no | yes (`atoms/sm.toml`, `ions/sm+2.toml`, `ions/sm+3.toml`) | done |  |
| Eu | yes | yes (+2,+3) | no | yes (`atoms/eu.toml`, `ions/eu+2.toml`, `ions/eu+3.toml`) | done |  |
| Gd | yes | yes (+3) | no | yes (`atoms/gd.toml`, `ions/gd+3.toml`) | done |  |
| Tb | yes | yes (+3,+4) | no | yes (`atoms/tb.toml`, `ions/tb+3.toml`, `ions/tb+4.toml`) | done |  |
| Dy | yes | yes (+3) | no | yes (`atoms/dy.toml`, `ions/dy+3.toml`) | done |  |
| Ho | yes | yes (+3) | no | yes (`atoms/ho.toml`, `ions/ho+3.toml`) | done |  |
| Er | yes | yes (+3) | no | yes (`atoms/er.toml`, `ions/er+3.toml`) | done |  |
| Tm | yes | yes (+3) | no | yes (`atoms/tm.toml`, `ions/tm+3.toml`) | done |  |
| Yb | yes | yes (+2,+3) | no | yes (`atoms/yb.toml`, `ions/yb+2.toml`, `ions/yb+3.toml`) | done |  |
| Lu | yes | yes (+3) | no | yes (`atoms/lu.toml`, `ions/lu+3.toml`) | done |  |
| Hf | yes | yes (+4) | no | yes (`atoms/hf.toml`, `ions/hf+4.toml`) | done |  |
| Ta | yes | no | no | yes (`atoms/ta.toml`) | done |  |
| W | yes | no | no | yes (`atoms/w.toml`) | done |  |
| Re | yes | no | no | yes (`atoms/re.toml`) | done |  |
| Os | yes | no | no | yes (`atoms/os.toml`) | done |  |
| Ir | yes | no | no | yes (`atoms/ir.toml`) | done |  |
| Pt | yes | yes (+2) | no | yes (`atoms/pt.toml`, `ions/pt+2.toml`) | done |  |
| Au | yes | yes (+1,+3) | no | yes (`atoms/au.toml`, `ions/au+1.toml`, `ions/au+3.toml`) | done | `Au3+` conformance pinned to low-spin |
| Hg | yes | yes (+2, [Hg2]2+) | no | yes (`atoms/hg.toml`, `ions/hg+2.toml`, `ions/hg+2_dimer.toml`) | done | Dimer encoded as two `Hg+1` atoms with single bond |
| Fr | yes | yes (+1) | no | yes (`atoms/fr.toml`, `ions/fr+1.toml`) | done |  |
| Ra | yes | yes (+2) | no | yes (`atoms/ra.toml`, `ions/ra+2.toml`) | done |  |
| Ac | yes | yes (+3) | no | yes (`atoms/ac.toml`, `ions/ac+3.toml`) | done |  |
| Th | yes | yes (+3,+4) | no | yes (`atoms/th.toml`, `ions/th+3.toml`, `ions/th+4.toml`) | done |  |
| Pa | yes | yes (+4,+5) | no | yes (`atoms/pa.toml`, `ions/pa+4.toml`, `ions/pa+5.toml`) | done | Valence table outer electrons adjusted to 5 |
| U | yes | yes (+3,+4,+5,+6) | no | yes (`atoms/u.toml`, `ions/u+3.toml`, `ions/u+4.toml`, `ions/u+5.toml`, `ions/u+6.toml`) | done | Valence table outer electrons adjusted to 6 |
| Np | yes | yes (+3,+4,+5,+6) | no | yes (`atoms/np.toml`, `ions/np+3.toml`, `ions/np+4.toml`, `ions/np+5.toml`, `ions/np+6.toml`) | done | Valence table outer electrons adjusted to 6 |
| Pu | yes | yes (+3,+4,+5,+6,+7) | no | yes (`atoms/pu.toml`, `ions/pu+3.toml`, `ions/pu+4.toml`, `ions/pu+5.toml`, `ions/pu+6.toml`, `ions/pu+7.toml`) | done | Valence table outer electrons adjusted to 7 |
| Am | yes | yes (+2,+3,+4) | no | yes (`atoms/am.toml`, `ions/am+2.toml`, `ions/am+3.toml`, `ions/am+4.toml`) | done |  |
| Cm | yes | yes (+3,+4) | no | yes (`atoms/cm.toml`, `ions/cm+3.toml`, `ions/cm+4.toml`) | done |  |
| Bk | yes | yes (+3,+4) | no | yes (`atoms/bk.toml`, `ions/bk+3.toml`, `ions/bk+4.toml`) | done |  |
| Cf | yes | yes (+2,+3,+4) | no | yes (`atoms/cf.toml`, `ions/cf+2.toml`, `ions/cf+3.toml`, `ions/cf+4.toml`) | done |  |
| Es | yes | yes (+2,+3) | no | yes (`atoms/es.toml`, `ions/es+2.toml`, `ions/es+3.toml`) | done |  |
| Fm | yes | yes (+2,+3) | no | yes (`atoms/fm.toml`, `ions/fm+2.toml`, `ions/fm+3.toml`) | done |  |
| Md | yes | yes (+2,+3) | no | yes (`atoms/md.toml`, `ions/md+2.toml`, `ions/md+3.toml`) | done |  |
| No | yes | yes (+2) | no | yes (`atoms/no.toml`, `ions/no+2.toml`) | done |  |
| Lr | yes | yes (+3) | no | yes (`atoms/lr.toml`, `ions/lr+3.toml`) | done |  |
| Rf | yes | no | no | yes (`atoms/rf.toml`) | done |  |
| Db | yes | no | no | yes (`atoms/db.toml`) | done |  |
| Sg | yes | no | no | yes (`atoms/sg.toml`) | done |  |
| Bh | yes | no | no | yes (`atoms/bh.toml`) | done |  |
| Hs | yes | no | no | yes (`atoms/hs.toml`) | done |  |

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
- `hydrides/`
- `covalent/`

## Current Conformance Set

- `atoms/`: `h`, `he`, `li`, `be`, `b`, `c`, `n`, `o`, `f`, `ne`, `na`, `mg`, `al`, `si`, `p`, `s`, `cl`, `ar`, `k`, `ca`, `sc`, `ti`, `v`, `cr`, `mn`, `fe`, `co`, `ni`, `cu`, `zn`, `ga`, `ge`, `as`, `se`, `br`, `kr`, `rb`, `sr`, `y`, `zr`, `nb`, `mo`, `tc`, `ru`, `rh`, `pd`, `ag`, `cd`, `in`, `sn`, `sb`, `te`, `i`, `xe`, `cs`, `ba`, `tl`, `pb`, `bi`, `po`, `at`, `rn`, `la`, `ce`, `pr`, `nd`, `pm`, `sm`, `eu`, `gd`, `tb`, `dy`, `ho`, `er`, `tm`, `yb`, `lu`, `hf`, `ta`, `w`, `re`, `os`, `ir`, `pt`, `au`, `hg`, `fr`, `ra`, `ac`, `th`, `pa`, `u`, `np`, `pu`, `am`, `cm`, `bk`, `cf`, `es`, `fm`, `md`, `no`, `lr`, `rf`, `db`, `sg`, `bh`, `hs`
- `ions/`: `h+1`, `li+1`, `be+2`, `mg+2`, `al+3`, `n-3`, `o-2`, `f-1`, `na+1`, `p-3`, `s-2`, `cl-1`, `k+1`, `ca+2`, `sc+3`, `ti+3`, `ti+4`, `v+2`, `v+3`, `cr+2`, `cr+3`, `mn+2`, `mn+3`, `fe+2`, `fe+3`, `co+2`, `co+3`, `ni+2`, `cu+1`, `cu+2`, `zn+2`, `ga+3`, `as-3`, `se-2`, `br-1`, `rb+1`, `sr+2`, `y+3`, `zr+4`, `nb+3`, `mo+2`, `mo+3`, `ru+2`, `ru+3`, `rh+2`, `rh+3`, `pd+2`, `ag+1`, `ag+2`, `cd+2`, `in+3`, `sb-3`, `te-2`, `i-1`, `cs+1`, `ba+2`, `tl+1`, `tl+3`, `pb+2`, `pb+4`, `bi+3`, `po-2`, `at-1`, `la+3`, `ce+3`, `ce+4`, `pr+3`, `pr+4`, `nd+3`, `pm+3`, `sm+2`, `sm+3`, `eu+2`, `eu+3`, `gd+3`, `tb+3`, `tb+4`, `dy+3`, `ho+3`, `er+3`, `tm+3`, `yb+2`, `yb+3`, `lu+3`, `hf+4`, `pt+2`, `au+1`, `au+3`, `hg+2`, `hg+2_dimer`, `fr+1`, `ra+2`, `ac+3`, `th+3`, `th+4`, `pa+4`, `pa+5`, `u+3`, `u+4`, `u+5`, `u+6`, `np+3`, `np+4`, `np+5`, `np+6`, `pu+3`, `pu+4`, `pu+5`, `pu+6`, `pu+7`, `am+2`, `am+3`, `am+4`, `cm+3`, `cm+4`, `bk+3`, `bk+4`, `cf+2`, `cf+3`, `cf+4`, `es+2`, `es+3`, `fm+2`, `fm+3`, `md+2`, `md+3`, `no+2`, `lr+3`
- `hydrides/`: `h2`, `bh3`, `ch4`, `nh3`, `h2o`, `hf`, `sih4`, `ph3`, `h2s`, `hcl`, `ash3`, `h2se`, `hbr`, `sbh3`, `h2te`, `hi`
- `covalent/`: `bf2`, `bcl4-`, `ch3`, `cf3-`, `ccl3+`, `ch2_singlet`, `ch2_triplet`, `sih3`, `sif3-`, `sicl3+`, `sih2_singlet`, `sih2_triplet`, `nh2`, `nh2-`, `nh4+`, `ph2`, `ph2-`, `ph4+`, `ash2`, `ash2-`, `ash4+`, `sbh2`, `sbh2-`, `sbh4+`, `oh`, `oh-`, `oh3+`, `sh`, `sh-`, `sh3+`, `seh`, `seh-`, `seh3+`, `teh`, `teh-`, `teh3+`, `f-`, `cl-`, `br-`, `i-`

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
- Added 6s/6p conformance coverage: atoms `Cs`, `Ba`, `Tl`, `Pb`, `Bi`, `Po`, `At`, `Rn`; ions `Cs+`, `Ba2+`, `Tl+`, `Tl3+`, `Pb2+`, `Pb4+`, `Bi3+`, `Po2-`, `At-`.
- Added 4f conformance coverage (`La`..`Lu`) with ions: all `3+`, plus `Ce4+`, `Pr4+`, `Tb4+`, and `Sm2+`, `Eu2+`, `Yb2+`.
- Added 5d conformance coverage (`Hf`..`Hg`) with ions: `Hf4+`, `Pt2+`, `Au+`, `Au3+`, `Hg2+`, and `[Hg2]2+` as `ions/hg+2_dimer.toml`.
- Pinned spin-sensitive transition-metal ion queries in conformance inputs (high-spin/low-spin where requested) and aligned `Pr` counts support (`default-valence-table`: `outer_electrons=5`, `allowed_valences=[3,4]`) so all conformance snapshots resolve with success.
- Added conformance coverage for `Fr`..`Lr` atoms and all corresponding registry ion charge states.
- Added neutral-atom conformance coverage for `Rf`, `Db`, `Sg`, `Bh`, `Hs`.
- Updated actinide counts input defaults by setting valence-table `outer_electrons`: `Pa=5`, `U=6`, `Np=6`, `Pu=7`.
- Added simple non-metal hydrides as atom-like conformance inputs (`?{EHn}` / `?{EH}`) with `implicit_h=true` and no explicit bonds.
- Added matching hydride atom-type specs to the default registry for atom-typing validation (`BH3`, `CH4`, `NH3`, `OH2`, `FH`, `SiH4`, `PH3`, `SH2`, `ClH`, `AsH3`, `SeH2`, `BrH`, `SbH3`, `TeH2`, `IH`).
- Counts-based typing now uses explicit bond-order sum for `v` (implicit hydrogens no longer inflate valence display), resolving outputs like `{O/2H2}` instead of `{O/2H2v2}`.
- Atoms/ions/hydrides rollout pass completed; table statuses are all `done`.
- Added covalent-state rollout phase in `covalent/` with `implicit_h=false` for all cases (all hydrogens explicit as atoms/bonds), following Hill+charge file naming (`CH3-.toml`, `NH4+.toml`) and `_isomer` suffix for spin-isomer cases (`CH2_singlet`, `CH2_triplet`).
- Covered required non-metal motif families with mixed substituents (`H`, `F`, `Cl`) across B/C/Si, pnictogens, chalcogens, and halogens, including radicals, anions, and cations.
- Added only required registry states for covalent motifs (`C/1v2`, `Si+v3`, `Si-/1v3`, `Si/1^2v2`, `P/1^1v2`, `P-/2v2`, `As/1^1v2`, `As+v4`, `As-/2v2`, `Sb/1^1v2`, `Sb+v4`, `Sb-/2v2`, `Se+/1v3`, `Te/2^1v1`, `Te+/1v3`, `Te-/3v1`).
- Conformance parity verified: `atom_typing` and `counts` both succeed on all new `covalent/` snapshots; no valence-table adjustments were required in this phase.
- Mixed implicit/explicit hydrogen representations are deferred to a future phase.