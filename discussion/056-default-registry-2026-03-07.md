# 056 — Default atom-type registry

Status: In Progress
Date: 2026-03-07
Relates: [217](217-rhea-participant-failures-2026-08-30.md)

The [HW extraction](#hw-extraction-spike-2026-10-31) records noble-gas evidence and
halogen candidates. The noble-gas audit covers earlier and newly added entries;
the supplied primary-paper PDFs have been reviewed. The noble-gas review and
selections are complete. Halogen implementation is paused while overlapping
representations are clarified.
The March planning and implementation notes below are historical; later updates to
the tracking table are dated.

This document tracks incremental, reviewable expansion of the default atom type registry.
Registry entries are treated as code: explicit scope, per-entry review, and matching conformance coverage.

## Project Status (2026-10-31)

**Last implementation check: 711 conformance cases and 366 valence unit tests passed.**
Each EDN input runs through AtomTyping and Counts with both Strict and MostSaturated policies. The selected
noble-gas HW additions and six xenon-anion replacements are implemented. Perxenate
uses Xe#v8 as an accepted assumption. All Rn compound entries, neutral Xe#m8
and the XeO4 adduct entries/fixtures are removed; MeCN→XeO3 uses only the Xe=O description.
Registry references are complete for this slice: HW, selected primary sources and PW assignment markers.
All 33 retained noble-gas states and the three F/N companion states have conformance coverage.
Overlapping halogen descriptions need clarification before implementation.

### Conformance suite breakdown

| Category | Files | Description |
|---|---|---|
| `atoms/` | 108 | Free atoms (ground/excited states), H through Hs |
| `atomic_ions/` | 127 | Atomic ions and Hg₂²⁺ |
| `hydrides/` | 17 | Hydrides with implicit H, plus explicit H₂ |
| `inorganic_small/` | 75 | Small inorganic molecules and salts |
| `inorganic_ions/` | 42 | Inorganic ions |
| `donor_acceptor/` | 4 | Donor–acceptor adducts with dative bonds |
| `hydrocarbons/` | 53 | C/H compounds |
| `functional_groups/` | 73 | Organic heteroatom groups (halides, amines, etc.) |
| `carbonyl_derivatives/` | 62 | Non-aromatic carbonyl compounds (aldehydes, ketones, acids, esters, amides, acyl derivatives) |
| `heterocycles/` | 53 | Non-aromatic O/N/S heterocycles, crown ethers |
| `aromatic/` | 34 | Aromatic carbocycles, heterocycles and ions |
| `hypervalent/` | 20 | Expanded-octet species (S, P, Cl, I) |
| `transition_metals/` | 7 | Transition-metal complexes |
| `stereo_tetrahedral/` | 22 | Tetrahedral stereochemistry |
| `stereo_cis_trans/` | 14 | Cis/trans stereochemistry |
| **Total** | **711** | |

## March 2026 record

### Key code files

- **Registry**: `umol-models-graph/config/default-registry.toml` — per-element atom type spec lists
- **Valence table**: `umol-models-graph/config/default-valence-table.toml` — element `outer_electrons` + `allowed_valences`
- **Test harness**: `umol-models-graph/tests/resolution/mod.rs` — rstest-based conformance runner with insta snapshots
- **Test data**: `umol-models-graph/tests/resolution/data/<category>/<compound>.toml`
- **Snapshots**: `umol-models-graph/tests/resolution/snapshots/`
- **Resolver**: `umol-models-graph/src/graph_ir/resolver.rs`
- **Valence logic**: `umol-models-graph/src/graph_ir/valence.rs`, `config_data.rs`
- **Atom type spec/query**: `umol-models-graph/src/graph_ir/atom_type.rs`
- **Tracking doc**: `discussion/56-default-registry-2026-03-07.md` (this file)

### Conformance execution notes

- Run the conformance suite with feature flag enabled:
  - `cargo test --features conformance --test resolution`
- Keep the compile-time refresh marker comment in `umol-models-graph/tests/resolution/mod.rs` up to date when adding/removing conformance files (current marker: `v6`).
- Validation criterion for accepted conformance snapshots:
  - `atom_typing.success: true`
  - `counts.success: true`

## Tracking Table

Noble-gas and F/N companion rows updated 2026-10-31; other rows retain the March record.

Status values:

- `todo`: not started
- `in_progress`: currently being edited
- `review`: implemented, awaiting chemistry review
- `done`: accepted

| Element | Ground state | Discrete ions | Typical valences | Conformance added | Status | Notes |
|---|---|---|---|---|---|---|
| H | yes | yes (+1) | yes | yes (`atoms/h.toml`, `atomic_ions/h+1.toml`, hydrides) | done |  |
| He | yes | no | n/a | yes (`atoms/he.edn`) | done | HW review complete; zero-lone-pair description retained |
| Li | yes | yes (+1) | no | yes (`atoms/li.toml`, `atomic_ions/li+1.toml`) | done |  |
| Be | yes | yes (+2) | no | yes (`atoms/be.toml`, `atomic_ions/be+2.toml`) | done |  |
| B | yes | no | no | yes (`atoms/b.toml`) | done |  |
| C | yes | no | no | yes (`atoms/c.toml`) | done |  |
| N | yes | yes (-3) | partial | yes (`atoms/n.edn`, `atomic_ions/n-3.edn`, `donor_acceptor/mecn-xeo3.edn`) | done | Noble-gas companion addition: N#v3#d |
| O | yes | yes (-2) | no | yes (`atoms/o.toml`, `atomic_ions/o-2.toml`) | done |  |
| F | yes | yes (-1) | partial | yes (`atoms/f.edn`, `atomic_ions/f-.edn`, hydrides, `inorganic_small/xef2.edn`, `donor_acceptor/xe2f11+.edn`) | done | Noble-gas companion additions: F#n3#m and F#c-#n2#d2; halogen review remains open |
| Ne | yes | no | n/a | yes (`atoms/ne.edn`) | done | HW review complete; zero-lone-pair description retained |
| Na | yes | yes (+1) | no | yes (`atoms/na.toml`, `atomic_ions/na+1.toml`) | done |  |
| Mg | yes | yes (+2) | no | yes (`atoms/mg.toml`, `atomic_ions/mg+2.toml`) | done |  |
| Al | yes | yes (+3) | no | yes (`atoms/al.toml`, `atomic_ions/al+3.toml`) | done |  |
| Si | yes | no | no | yes (`atoms/si.toml`) | done |  |
| P | yes | yes (-3) | no | yes (`atoms/p.toml`, `atomic_ions/p-3.toml`) | done |  |
| S | yes | yes (-2) | no | yes (`atoms/s.toml`, `atomic_ions/s-2.toml`) | done |  |
| Cl | yes | yes (-1) | no | yes (`atoms/cl.toml`, `atomic_ions/cl-1.toml`) | done |  |
| Ar | yes | no | n/a | yes (`atoms/ar.edn`) | done | HW review complete; zero-lone-pair description retained |
| K | yes | yes (+1) | no | yes (`atoms/k.toml`, `atomic_ions/k+1.toml`) | done |  |
| Ca | yes | yes (+2) | no | yes (`atoms/ca.toml`, `atomic_ions/ca+2.toml`) | done |  |
| Sc | yes | yes (+3) | no | yes (`atoms/sc.toml`, `atomic_ions/sc+3.toml`) | done |  |
| Ti | yes | yes (+3,+4) | no | yes (`atoms/ti.toml`, `atomic_ions/ti+3.toml`, `atomic_ions/ti+4.toml`) | done |  |
| V | yes | yes (+2,+3) | no | yes (`atoms/v.toml`, `atomic_ions/v+2.toml`, `atomic_ions/v+3.toml`) | done |  |
| Cr | yes | yes (+2,+3) | no | yes (`atoms/cr.toml`, `atomic_ions/cr+2.toml`, `atomic_ions/cr+3.toml`) | done |  |
| Mn | yes | yes (+2,+3) | no | yes (`atoms/mn.toml`, `atomic_ions/mn+2.toml`, `atomic_ions/mn+3.toml`) | done |  |
| Fe | yes | yes (+2,+3) | no | yes (`atoms/fe.toml`, `atomic_ions/fe+2.toml`, `atomic_ions/fe+3.toml`) | done |  |
| Co | yes | yes (+2,+3) | no | yes (`atoms/co.toml`, `atomic_ions/co+2.toml`, `atomic_ions/co+3.toml`) | done |  |
| Ni | yes | yes (+2) | no | yes (`atoms/ni.toml`, `atomic_ions/ni+2.toml`) | done |  |
| Cu | yes | yes (+1,+2) | no | yes (`atoms/cu.toml`, `atomic_ions/cu+1.toml`, `atomic_ions/cu+2.toml`) | done |  |
| Zn | yes | yes (+2) | no | yes (`atoms/zn.toml`, `atomic_ions/zn+2.toml`) | done |  |
| Ga | yes | yes (+3) | no | yes (`atoms/ga.toml`, `atomic_ions/ga+3.toml`) | done |  |
| Ge | yes | no | no | yes (`atoms/ge.toml`) | done |  |
| As | yes | yes (-3) | no | yes (`atoms/as.toml`, `atomic_ions/as-3.toml`) | done |  |
| Se | yes | yes (-2) | no | yes (`atoms/se.toml`, `atomic_ions/se-2.toml`) | done |  |
| Br | yes | yes (-1) | no | yes (`atoms/br.toml`, `atomic_ions/br-1.toml`, hydrides) | done |  |
| Kr | yes | no | yes (selected) | yes (`atoms/kr.edn`, `inorganic_small/krf2.edn`, `inorganic_ions/krf+.edn`, `donor_acceptor/hcn-krf+.edn`) | done | HW review complete; all four retained states covered |
| Rb | yes | yes (+1) | no | yes (`atoms/rb.toml`, `atomic_ions/rb+1.toml`) | done |  |
| Sr | yes | yes (+2) | no | yes (`atoms/sr.toml`, `atomic_ions/sr+2.toml`) | done |  |
| Y | yes | yes (+3) | no | yes (`atoms/y.toml`, `atomic_ions/y+3.toml`) | done |  |
| Zr | yes | yes (+4) | no | yes (`atoms/zr.toml`, `atomic_ions/zr+4.toml`) | done |  |
| Nb | yes | yes (+3) | no | yes (`atoms/nb.toml`, `atomic_ions/nb+3.toml`) | done |  |
| Mo | yes | yes (+2,+3) | no | yes (`atoms/mo.toml`, `atomic_ions/mo+2.toml`, `atomic_ions/mo+3.toml`) | done |  |
| Tc | yes | no | no | yes (`atoms/tc.toml`) | done | No ion added in current scope |
| Ru | yes | yes (+2,+3) | no | yes (`atoms/ru.toml`, `atomic_ions/ru+2.toml`, `atomic_ions/ru+3.toml`) | done |  |
| Rh | yes | yes (+2,+3) | no | yes (`atoms/rh.toml`, `atomic_ions/rh+2.toml`, `atomic_ions/rh+3.toml`) | done |  |
| Pd | yes | yes (+2) | no | yes (`atoms/pd.toml`, `atomic_ions/pd+2.toml`) | done |  |
| Ag | yes | yes (+1,+2) | no | yes (`atoms/ag.toml`, `atomic_ions/ag+1.toml`, `atomic_ions/ag+2.toml`) | done |  |
| Cd | yes | yes (+2) | no | yes (`atoms/cd.toml`, `atomic_ions/cd+2.toml`) | done |  |
| In | yes | yes (+3) | no | yes (`atoms/in.toml`, `atomic_ions/in+3.toml`) | done |  |
| Sn | yes | no | no | yes (`atoms/sn.toml`) | done |  |
| Sb | yes | yes (-3) | no | yes (`atoms/sb.toml`, `atomic_ions/sb-3.toml`) | done |  |
| Te | yes | yes (-2) | no | yes (`atoms/te.toml`, `atomic_ions/te-2.toml`) | done |  |
| I | yes | yes (-1) | no | yes (`atoms/i.toml`, `atomic_ions/i-1.toml`, hydrides) | done |   |
| Xe | yes | no | yes (selected) | yes (`atoms/xe.edn`, `inorganic_small/xe*.edn`, `inorganic_ions/xe*.edn`, `donor_acceptor/`) | done | HW review complete; all 25 retained states covered; HW/primary/PW references added |
| Cs | yes | yes (+1) | no | yes (`atoms/cs.toml`, `atomic_ions/cs+1.toml`) | done |  |
| Ba | yes | yes (+2) | no | yes (`atoms/ba.toml`, `atomic_ions/ba+2.toml`) | done |  |
| Tl | yes | yes (+1,+3) | no | yes (`atoms/tl.toml`, `atomic_ions/tl+1.toml`, `atomic_ions/tl+3.toml`) | done |  |
| Pb | yes | yes (+2,+4) | no | yes (`atoms/pb.toml`, `atomic_ions/pb+2.toml`, `atomic_ions/pb+4.toml`) | done |  |
| Bi | yes | yes (+3) | no | yes (`atoms/bi.toml`, `atomic_ions/bi+3.toml`) | done |  |
| Po | yes | yes (-2) | no | yes (`atoms/po.toml`, `atomic_ions/po-2.toml`) | done |  |
| At | yes | yes (-1) | no | yes (`atoms/at.toml`, `atomic_ions/at-1.toml`) | done |  |
| Rn | yes | no | none selected | yes (`atoms/rn.edn`) | done | HW review complete; unsupported compound states removed; Rn#n4 retained |
| La | yes | yes (+3) | no | yes (`atoms/la.toml`, `atomic_ions/la+3.toml`) | done |  |
| Ce | yes | yes (+3,+4) | no | yes (`atoms/ce.toml`, `atomic_ions/ce+3.toml`, `atomic_ions/ce+4.toml`) | done |  |
| Pr | yes | yes (+3,+4) | no | yes (`atoms/pr.toml`, `atomic_ions/pr+3.toml`, `atomic_ions/pr+4.toml`) | done | Counts requires `Pr` valence-table support for +4 |
| Nd | yes | yes (+3) | no | yes (`atoms/nd.toml`, `atomic_ions/nd+3.toml`) | done |  |
| Pm | yes | yes (+3) | no | yes (`atoms/pm.toml`, `atomic_ions/pm+3.toml`) | done |  |
| Sm | yes | yes (+2,+3) | no | yes (`atoms/sm.toml`, `atomic_ions/sm+2.toml`, `atomic_ions/sm+3.toml`) | done |  |
| Eu | yes | yes (+2,+3) | no | yes (`atoms/eu.toml`, `atomic_ions/eu+2.toml`, `atomic_ions/eu+3.toml`) | done |  |
| Gd | yes | yes (+3) | no | yes (`atoms/gd.toml`, `atomic_ions/gd+3.toml`) | done |  |
| Tb | yes | yes (+3,+4) | no | yes (`atoms/tb.toml`, `atomic_ions/tb+3.toml`, `atomic_ions/tb+4.toml`) | done |  |
| Dy | yes | yes (+3) | no | yes (`atoms/dy.toml`, `atomic_ions/dy+3.toml`) | done |  |
| Ho | yes | yes (+3) | no | yes (`atoms/ho.toml`, `atomic_ions/ho+3.toml`) | done |  |
| Er | yes | yes (+3) | no | yes (`atoms/er.toml`, `atomic_ions/er+3.toml`) | done |  |
| Tm | yes | yes (+3) | no | yes (`atoms/tm.toml`, `atomic_ions/tm+3.toml`) | done |  |
| Yb | yes | yes (+2,+3) | no | yes (`atoms/yb.toml`, `atomic_ions/yb+2.toml`, `atomic_ions/yb+3.toml`) | done |  |
| Lu | yes | yes (+3) | no | yes (`atoms/lu.toml`, `atomic_ions/lu+3.toml`) | done |  |
| Hf | yes | yes (+4) | no | yes (`atoms/hf.toml`, `atomic_ions/hf+4.toml`) | done |  |
| Ta | yes | no | no | yes (`atoms/ta.toml`) | done |  |
| W | yes | no | no | yes (`atoms/w.toml`) | done |  |
| Re | yes | no | no | yes (`atoms/re.toml`) | done |  |
| Os | yes | no | no | yes (`atoms/os.toml`) | done |  |
| Ir | yes | no | no | yes (`atoms/ir.toml`) | done |  |
| Pt | yes | yes (+2) | no | yes (`atoms/pt.toml`, `atomic_ions/pt+2.toml`) | done |  |
| Au | yes | yes (+1,+3) | no | yes (`atoms/au.toml`, `atomic_ions/au+1.toml`, `atomic_ions/au+3.toml`) | done | `Au3+` conformance pinned to low-spin |
| Hg | yes | yes (+2, [Hg2]2+) | no | yes (`atoms/hg.toml`, `atomic_ions/hg+2.toml`, `atomic_ions/hg+2_dimer.toml`) | done | Dimer encoded as two `Hg+1` atoms with single bond |
| Fr | yes | yes (+1) | no | yes (`atoms/fr.toml`, `atomic_ions/fr+1.toml`) | done |  |
| Ra | yes | yes (+2) | no | yes (`atoms/ra.toml`, `atomic_ions/ra+2.toml`) | done |  |
| Ac | yes | yes (+3) | no | yes (`atoms/ac.toml`, `atomic_ions/ac+3.toml`) | done |  |
| Th | yes | yes (+3,+4) | no | yes (`atoms/th.toml`, `atomic_ions/th+3.toml`, `atomic_ions/th+4.toml`) | done |  |
| Pa | yes | yes (+4,+5) | no | yes (`atoms/pa.toml`, `atomic_ions/pa+4.toml`, `atomic_ions/pa+5.toml`) | done | Valence table outer electrons adjusted to 5 |
| U | yes | yes (+3,+4,+5,+6) | no | yes (`atoms/u.toml`, `atomic_ions/u+3.toml`, `atomic_ions/u+4.toml`, `atomic_ions/u+5.toml`, `atomic_ions/u+6.toml`) | done | Valence table outer electrons adjusted to 6 |
| Np | yes | yes (+3,+4,+5,+6) | no | yes (`atoms/np.toml`, `atomic_ions/np+3.toml`, `atomic_ions/np+4.toml`, `atomic_ions/np+5.toml`, `atomic_ions/np+6.toml`) | done | Valence table outer electrons adjusted to 6 |
| Pu | yes | yes (+3,+4,+5,+6,+7) | no | yes (`atoms/pu.toml`, `atomic_ions/pu+3.toml`, `atomic_ions/pu+4.toml`, `atomic_ions/pu+5.toml`, `atomic_ions/pu+6.toml`, `atomic_ions/pu+7.toml`) | done | Valence table outer electrons adjusted to 7 |
| Am | yes | yes (+2,+3,+4) | no | yes (`atoms/am.toml`, `atomic_ions/am+2.toml`, `atomic_ions/am+3.toml`, `atomic_ions/am+4.toml`) | done |  |
| Cm | yes | yes (+3,+4) | no | yes (`atoms/cm.toml`, `atomic_ions/cm+3.toml`, `atomic_ions/cm+4.toml`) | done |  |
| Bk | yes | yes (+3,+4) | no | yes (`atoms/bk.toml`, `atomic_ions/bk+3.toml`, `atomic_ions/bk+4.toml`) | done |  |
| Cf | yes | yes (+2,+3,+4) | no | yes (`atoms/cf.toml`, `atomic_ions/cf+2.toml`, `atomic_ions/cf+3.toml`, `atomic_ions/cf+4.toml`) | done |  |
| Es | yes | yes (+2,+3) | no | yes (`atoms/es.toml`, `atomic_ions/es+2.toml`, `atomic_ions/es+3.toml`) | done |  |
| Fm | yes | yes (+2,+3) | no | yes (`atoms/fm.toml`, `atomic_ions/fm+2.toml`, `atomic_ions/fm+3.toml`) | done |  |
| Md | yes | yes (+2,+3) | no | yes (`atoms/md.toml`, `atomic_ions/md+2.toml`, `atomic_ions/md+3.toml`) | done |  |
| No | yes | yes (+2) | no | yes (`atoms/no.toml`, `atomic_ions/no+2.toml`) | done |  |
| Lr | yes | yes (+3) | no | yes (`atoms/lr.toml`, `atomic_ions/lr+3.toml`) | done |  |
| Rf | yes | no | no | yes (`atoms/rf.toml`) | done |  |
| Db | yes | no | no | yes (`atoms/db.toml`) | done |  |
| Sg | yes | no | no | yes (`atoms/sg.toml`) | done |  |
| Bh | yes | no | no | yes (`atoms/bh.toml`) | done |  |
| Hs | yes | no | no | yes (`atoms/hs.toml`) | done |  |

Extend this table as new elements are considered.

### Organic compounds

| Class | Scope | Conformance dir | Status | Notes |
|-------|-------|-----------------|--------|-------|
| Hydrocarbons | Alkanes, alkenes, alkynes, dienes/cumulenes, cycloalkanes, unsaturated alicycles, fused/spiro/bridged bicycles, polycycles | `hydrocarbons/` | done | 51 files, `implicit_h=true`, `?{CHn}` |
| Alcohols / ethers | R–OH, R–O–R′ | — | todo | |
| Carbonyl derivatives | Aldehydes, ketones, carboxylic acids, esters, amides | `carbonyl_derivatives/` | review | 62 files, `implicit_h=true`, non-aromatic |
| Organic halides | R–F, R–Cl, R–Br, R–I | — | todo | |
| Other | Amines, thiols, etc. | — | todo | |
| Functional groups | Halides, azides, cyanates, thiocyanates, alcohols, ethers, amines, nitroso, nitro, sulfoxides, hydrazines, diazo/azo | `functional_groups/` | review | 73 files, `implicit_h=true` |
| Heterocycles | Non-aromatic O/N/S heterocycles, crown ethers, bridged/fused | `heterocycles/` | review | 52 files, `implicit_h=true` |

Extend this table as new organic classes are added.

## Conformance Policy

Per element, add explicit resolution inputs that cover:

- Free atom query (with explicit `implicit_h` policy)
- At least one bonded case when chemically meaningful
- Relevant charge states (for ion stage)
- At least one "should fail" case where useful to constrain over-broad matching

Suggested naming pattern:

- Atomic ground state: `tests/resolution/data/atoms/<element>.toml`
- Atomic valence/excited variants: `tests/resolution/data/atoms/<element>_<motif>.toml`
- Ionic ground state: `tests/resolution/data/atomic_ions/<element><charge>.toml`
- Ionic valence/excited variants: `tests/resolution/data/atomic_ions/<element><charge>_<motif>.toml`

Conventions:

- `<element>`: lowercase element symbol/name token used in file naming.
- `<charge>`: signed integer with explicit sign, e.g. `+1`, `-2` (no underscore between element and charge).
- `<motif>`: lowercase, underscore-separated descriptor; may encode occupation and optional multiplet.

## Categories

The [conformance suite breakdown](#conformance-suite-breakdown) lists the current
categories. Donor–acceptor adducts use `donor_acceptor/`, including charged adducts.

The former `implicit_h=true/false` setting was the precursor to atom defaults.
EDN inputs use MoleculeDefaults with per-input `:config-overrides`; hydrogen
handling follows those defaults and explicit atom predicates.

## Current Conformance Set

- `atoms/`: `h`, `he`, `li`, `be`, `b`, `c`, `n`, `o`, `f`, `ne`, `na`, `mg`, `al`, `si`, `p`, `s`, `cl`, `ar`, `k`, `ca`, `sc`, `ti`, `v`, `cr`, `mn`, `fe`, `co`, `ni`, `cu`, `zn`, `ga`, `ge`, `as`, `se`, `br`, `kr`, `rb`, `sr`, `y`, `zr`, `nb`, `mo`, `tc`, `ru`, `rh`, `pd`, `ag`, `cd`, `in`, `sn`, `sb`, `te`, `i`, `xe`, `cs`, `ba`, `tl`, `pb`, `bi`, `po`, `at`, `rn`, `la`, `ce`, `pr`, `nd`, `pm`, `sm`, `eu`, `gd`, `tb`, `dy`, `ho`, `er`, `tm`, `yb`, `lu`, `hf`, `ta`, `w`, `re`, `os`, `ir`, `pt`, `au`, `hg`, `fr`, `ra`, `ac`, `th`, `pa`, `u`, `np`, `pu`, `am`, `cm`, `bk`, `cf`, `es`, `fm`, `md`, `no`, `lr`, `rf`, `db`, `sg`, `bh`, `hs`
- `atomic_ions/`: `h+1`, `li+1`, `be+2`, `mg+2`, `al+3`, `n-3`, `o-2`, `f-1`, `na+1`, `p-3`, `s-2`, `cl-1`, `k+1`, `ca+2`, `sc+3`, `ti+3`, `ti+4`, `v+2`, `v+3`, `cr+2`, `cr+3`, `mn+2`, `mn+3`, `fe+2`, `fe+3`, `co+2`, `co+3`, `ni+2`, `cu+1`, `cu+2`, `zn+2`, `ga+3`, `as-3`, `se-2`, `br-1`, `rb+1`, `sr+2`, `y+3`, `zr+4`, `nb+3`, `mo+2`, `mo+3`, `ru+2`, `ru+3`, `rh+2`, `rh+3`, `pd+2`, `ag+1`, `ag+2`, `cd+2`, `in+3`, `sb-3`, `te-2`, `i-1`, `cs+1`, `ba+2`, `tl+1`, `tl+3`, `pb+2`, `pb+4`, `bi+3`, `po-2`, `at-1`, `la+3`, `ce+3`, `ce+4`, `pr+3`, `pr+4`, `nd+3`, `pm+3`, `sm+2`, `sm+3`, `eu+2`, `eu+3`, `gd+3`, `tb+3`, `tb+4`, `dy+3`, `ho+3`, `er+3`, `tm+3`, `yb+2`, `yb+3`, `lu+3`, `hf+4`, `pt+2`, `au+1`, `au+3`, `hg+2`, `hg+2_dimer`, `fr+1`, `ra+2`, `ac+3`, `th+3`, `th+4`, `pa+4`, `pa+5`, `u+3`, `u+4`, `u+5`, `u+6`, `np+3`, `np+4`, `np+5`, `np+6`, `pu+3`, `pu+4`, `pu+5`, `pu+6`, `pu+7`, `am+2`, `am+3`, `am+4`, `cm+3`, `cm+4`, `bk+3`, `bk+4`, `cf+2`, `cf+3`, `cf+4`, `es+2`, `es+3`, `fm+2`, `fm+3`, `md+2`, `md+3`, `no+2`, `lr+3`
- `hydrides/`: `h2`, `bh3`, `ch4`, `nh3`, `h2o`, `hf`, `sih4`, `ph3`, `h2s`, `hcl`, `ash3`, `h2se`, `hbr`, `sbh3`, `h2te`, `hi`
- `inorganic_small/`: neutral small inorganics (explicit H only): e.g. `bf2`, `bf3`, `bn`, `c2`, `ch3`, `ch2_singlet`, `ch2_triplet`, `sih3`, `sih2_*`, `nh2`, `ph2`, `ash2`, `sbh2`, `oh`, `sh`, `seh`, `teh`, `br2`, `cl2`, `cn-`, `f2`, `i2`, `n2`, `o2`, etc. Added nitrogen oxides (NO, NO2, N2O, N2O3, N2O4), carbon compounds (CO2, CS2, COS, HCN), oxyacids/hydroxides (HNO2, HNO3, HClO, H2CO3, H3BO3, HOBr, HOI), N-H compounds (N2H4, NH2OH, N2H2), peroxides (H2O2, H2S2), halides (NF3, NCl3, PF3, PCl3, SF2, SiF4, SiCl4, BCl3, ClF, BrF), and other (O3, HNCO, HOCN, HSCN, HNCS).
- `inorganic_ions/`: molecular ions (explicit H only): e.g. `bcl4-`, `cf3-`, `ccl3+`, `nh2-`, `nh4+`, `ph2-`, `ph4+`, `ash2-`, `ash4+`, `sbh2-`, `sbh4+`, `oh-`, `oh3+`, `sh-`, `sh3+`, `seh-`, `seh3+`, `teh-`, `teh3+`, `f-`, `cl-`, `br-`, `i-`, `sif3-`, `sicl3+`, etc.
- `hydrocarbons/`: `ethane`, `propane`, `butane`, `isobutane`, `pentane`, `isopentane`, `neopentane`, `hexane`, `ethene`, `propene`, `but-1-ene`, `but-2-ene`, `2-methylpropene`, `2-methylbut-2-ene`, `ethyne`, `propyne`, `but-1-yne`, `but-2-yne`, `allene`, `buta-1,3-diene`, `penta-1,4-diene`, `butatriene`, `cyclopropane`, `cyclobutane`, `cyclopentane`, `cyclohexane`, `cycloheptane`, `cyclopropene`, `cyclobutene`, `cyclopentene`, `cyclopentadiene`, `cyclohexene`, `cyclohexa-1,3-diene`, `cyclohexa-1,4-diene`, `bicyclo-1.1.0-butane`, `bicyclo-2.1.0-pentane`, `bicyclo-2.2.0-hexane`, `bicyclo-3.3.0-octane`, `bicyclo-4.3.0-nonane`, `decalin`, `spiropentane`, `spirohexane`, `spiroheptane`, `spirononane`, `bicyclo-1.1.1-pentane`, `norbornane`, `bicyclo-2.2.2-octane`, `tetrahedrane`, `prismane`, `cubane`, `adamantane`
- `functional_groups/`: methyl/ethyl × {fluoride, chloride, bromide, iodide, azide, isocyanide, cyanate, thiocyanate, isothiocyanate, nitrite, nitrate, thiol, selenol} (26); {propyl, isopropyl, butyl, isobutyl, tert-butyl, cyclopropyl, cyclohexyl}-chloride (7); dichloromethane, chloroform, carbon-tetrachloride, 1,1-dichloroethane, 1,2-dichloroethane, 1,2-dichlorocyclohexane, 1-chloro-2-hydroxyethane (7); methanol, ethanol, propan-1-ol, propan-2-ol, butan-1-ol, butan-2-ol, 2,2-dimethylpropan-2-ol, hydroxycyclopropane, hydroxycyclohexane (9); dimethylether, diethylether, 1,2-dimethoxyethane (3); methylamine, ethylamine, dimethylamine, trimethylamine, tetramethylammonium, trimethylamine-oxide (6); nitrosomethane, nitrosoethane, nitromethane, nitroethane, nitrosocyclohexane, nitrocyclohexane (6); dimethyl-sulfoxide, methylhydrazine, 1,1-dimethylhydrazine, 1,2-dimethylhydrazine, diazomethane, dimethyldiazene (6). Total: 73 files, `implicit_h=true`.
- `carbonyl_derivatives/`: non-aromatic carbonyl compounds (all `implicit_h=true`): aldehydes, ketones, dicarbonyls/hydroxy-carbonyls, monocarboxylic and dicarboxylic acids, acyl halides, anhydrides, carboxylates, esters/lactones/thioesters, amides/lactams/imides, and carbonyl-derived C=N species (oxime/hydrazone), plus pyruvic acid. Total: 62 files.
- `hypervalent/`: expanded-octet species (all `implicit_h=false`): SO2, SO3, H2SO3, H2SO4, SF4, SF6, H3PO3, H3PO4, P4O10, PF5, PCl5, POCl3, HClO2, HClO3, HClO4, ClF3, IF5, SO4(2-), PO4(3-), ClO4(-). Total: 20 files.
- `heterocycles/`: non-aromatic heterocycles (all `implicit_h=true`). O-saturated (oxirane, dioxirane, oxetane, 1,2-dioxetane, 1,3-dioxetane, tetrahydrofuran, tetrahydropyran); N-saturated (aziridine, azetidine, pyrrolidine, piperidine); S-saturated (thiirane, thietane, tetrahydrothiophene, thiane); O/S-unsaturated (2H-oxete, 2H-thiete, 2H-pyran, 4H-pyran, 2H-thiopyran, 4H-thiopyran); N-unsaturated (2H-azirine, 2,3-dihydroazete, azete, 3-pyrroline, 2-pyrroline); functionalized (sulfolane, 2-oxazolidone, hydantoin); multi-O/S (1,3-dioxolane, 1,2,3-trioxolane, 1,4-dioxane, 1,3-dithiane, 1,4-dithiane, 1,3,5-trithiane); multi-N (1,3-diazetidine, pyrazolidine, imidazolidine, 2-pyrazoline, 2-imidazoline); mixed (1,2-oxazetidine, 1,2-oxathiolane, 1,3-oxathiolane, morpholine, thiomorpholine); bridged/fused (pyrrolizidine, quinuclidine, 1-azaadamantane); crown ethers (12-crown-4, 15-crown-5, 18-crown-6). Total: 52 files.

## Decision Log

### 2026-03-07

- Agreed staged rollout: ground states -> discrete ions -> typical non-metal valences.
- Agreed registry expansion must be paired with conformance additions.
- Agreed registry specs are reviewed with code-level rigor (small, explicit, traceable changes).
- Agreed conformance naming/categories: `atoms/`, `atomic_ions/`, `hydrides/` and file naming rules above.
- Added `implicit_h` control to conformance inputs; counts strategy now honors this via `enable_implicit_hydrogens`.
- Expanded coverage to period-3 atoms (`Na`..`Ar`) and matching ions (`Na+`, `Al3+`, `P3-`, `S2-`, `Cl-`), plus `Be2+`.
- Expanded coverage to period-4 and period-5 s/p blocks, 3d block, and 4d block with ions: `Y3+`, `Zr4+`, `Nb3+`, `Mo2+`, `Mo3+`, `Ru2+`, `Ru3+`, `Rh2+`, `Rh3+`, `Pd2+`, `Ag+`, `Ag2+`, `Cd2+`.
- Added missing 3d atom conformance files (`atoms/sc.toml` through `atoms/zn.toml`) and refreshed snapshots.
- Added 6s/6p conformance coverage: atoms `Cs`, `Ba`, `Tl`, `Pb`, `Bi`, `Po`, `At`, `Rn`; ions `Cs+`, `Ba2+`, `Tl+`, `Tl3+`, `Pb2+`, `Pb4+`, `Bi3+`, `Po2-`, `At-`.
- Added 4f conformance coverage (`La`..`Lu`) with ions: all `3+`, plus `Ce4+`, `Pr4+`, `Tb4+`, and `Sm2+`, `Eu2+`, `Yb2+`.
- Added 5d conformance coverage (`Hf`..`Hg`) with ions: `Hf4+`, `Pt2+`, `Au+`, `Au3+`, `Hg2+`, and `[Hg2]2+` as `atomic_ions/hg+2_dimer.toml`.
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
- Added `hydrocarbons/` conformance category (51 files) with `implicit_h=true` using `?{CHn}` atom queries and C-C bonds only. Covers alkanes (8), alkenes (6), alkynes (4), dienes/cumulenes (4), cycloalkanes (5), unsaturated alicycles (7), fused bicycles (6), spiro bicycles (4), bridged bicycles (3), and polycycles (4: tetrahedrane, prismane, cubane, adamantane). No new registry entries required. All 51 tests pass for both `atom_typing` and `counts` strategies.
- Restructured conformance categories: renamed `ions/` → `atomic_ions/` (single-atom ions). Split `covalent/` → `inorganic_small/` (neutral small inorganics, including sextet species) and `inorganic_ions/` (molecular ions). Policy: `hydrides/` use implicit H only (one-atom graphs); `inorganic_small/`, `inorganic_ions/`, and hydrocarbons use explicit H where present (hydrocarbons keep implicit H in queries). Compromise between fine-grained compound classes and a single mixed category.
- Added `functional_groups/` conformance category (64 files) with `implicit_h=true`. Covers organic halides (methyl/ethyl × 4 halogens + 7 alkyl chlorides + 7 poly/mixed chloro compounds = 39), azides, isocyanides, cyanates, thiocyanates, isothiocyanates, nitrites, nitrates, thiols, selenols (methyl/ethyl × 9 = 18), alcohols (9), ethers (3), amines (6), nitroso (3), nitro (3). Inputs reviewed: formal charges on N+/O-/C-/N- where required (nitro, nitrate, azide, isocyanide, TMAO, tetramethylammonium). All 64 counts-strategy results pass. 26 atom-typing results show `ValenceAmbiguous` (registry disambiguation work remaining).
- Fixed `compute_implicit_hydrogens` to use RDKit-style effective-element lookup via `element.shift(-charge)` for allowed-valence comparison. Previously, charged atoms with specific allowed valences (e.g. N+ needing C's `[4]` instead of N's `[3]`) would fail to match.
- Extended Cl `allowed_valences` to `[1, 3, 5, 7]` in `default-valence-table.toml` to support hypervalent chlorine compounds (HClO2, HClO3, HClO4, ClF3).
- Added 36 new `inorganic_small/` conformance files: nitrogen oxides (NO, NO2, N2O, N2O3, N2O4), carbon compounds (CO2, CS2, COS, HCN), oxyacids/hydroxides (HNO2, HNO3, HClO, H2CO3, H3BO3, HOBr, HOI), N-H compounds (N2H4, NH2OH, N2H2), peroxides (H2O2, H2S2), halides (NF3, NCl3, PF3, PCl3, SF2, SiF4, SiCl4, BCl3, ClF, BrF), other (O3, HNCO, HOCN, HSCN, HNCS). All `implicit_h=false` with explicit H atoms and formal charges.
- Added 17 new `hypervalent/` conformance files: sulfur (SO2, SO3, H2SO3, H2SO4, SF4, SF6), phosphorus (H3PO3, H3PO4, P4O10, PF5, PCl5, POCl3), chlorine (HClO2, HClO3, HClO4, ClF3), iodine (IF5). All `implicit_h=false`.
- Resolution results: all 53 new files resolve successfully with `counts` strategy. `atom_typing` fails on 3 files (NO2, IF5, HClO4) due to missing registry entries for radical N (`{N/1^1v3}`), hypervalent I (`{I/1v5}`), and heptavalent Cl (`{Clv7}`). N2O4 input was corrected to symmetric O=N+(O-)–N+(O-)=O structure; both strategies now succeed.
- Identified 5 new atom type specs needed for full `atom_typing` coverage of new compounds: `{Cl/2v3}` (Cl valence 3), `{Cl/1v5}` (Cl valence 5), `{Clv7}` (Cl valence 7), `{I/1v5}` (I valence 5), `{N/1^1v3}` (radical N in NO2).
- Added `heterocycles/` conformance category (53 files) with `implicit_h=true`. Covers O-saturated (7), N-saturated (5), S-saturated (4), O/S-unsaturated (6), N-unsaturated (5), functionalized (4: sulfolane, succinimide, 2-oxazolidone, hydantoin), multi-O/S (6), multi-N (5), mixed heteroatom (5), bridged/fused (3: pyrrolizidine, quinuclidine, 1-azaadamantane), crown ethers (3: 12-crown-4, 15-crown-5, 18-crown-6). All 53 counts-strategy results pass. 5 atom-typing results show `ValenceAmbiguous` on carbonyl O (same registry disambiguation issue as functional_groups).
- Added new `carbonyl_derivatives/` conformance category (62 files, non-aromatic) covering aldehydes, ketones, acids, esters, amides, acyl halides, anhydrides, carboxylates, lactones/thioesters, and C=N derivatives tied to carbonyl chemistry.
- Added 6 `functional_groups/` compounds requested for S/N chemistry expansion: dimethyl-sulfoxide, methylhydrazine, 1,1-dimethylhydrazine, 1,2-dimethylhydrazine, diazomethane, dimethyldiazene.
- Moved `beta-lactam` and `succinimide` from `heterocycles/` into `carbonyl_derivatives/` to keep carbonyl-centered compounds in one category.

## HW extraction spike (2026-10-31)

For each description, distinguish evidence for the species from evidence for the
bonding used to construct its encoding. Where HW is insufficient, follow its
references to primary papers (through reviews where necessary), or search independently.
Record the supporting passage and flag unresolved assignments. Electron conservation
and conformance checks do not establish source support for a bonding description.

Registry comments use HW volume/page references as the baseline. A source table at
the top of the registry defines short keys; primary references supplement HW only
where they add information used in an assignment. PW means “Present work —
calculations and adopted assignments,” marked as `assignment: PW` on the relevant
examples. Methods, results and assumptions are recorded here.

### Noble Gases

Source: [Holleman/Wiberg 103, Band 1, chapter XI](<../materials/valence/Holleman-Wiberg - Anorganische Chemie - 103. Auflage - Band 1.pdf#page=497>),
printed pp. 459–474, PDF pp. 497–512. Page references below are printed pages.
Atom strings use AtomTypeRegistry::raise_defaults. This audit covers all **33 current
noble-gas entries**, plus the F/N companions. Earlier and newly added entries use the
same evidence criterion. These are evidence findings; registry selections remain with
the user.

Current tally: **27 supported**, **6 inferred** and **0 open**.
F#n3#m and N#v3#d are supported. F#c-#n2#d2 is also selected for the
Xe2F11+ bridge as a constructed donor–acceptor description.

#### Current entries: evidence review

**Supported** means the indicated source description supports the constructed atom
state. **Inferred** means selected under the construction convention below, without
a verified individual bonding assignment. **Open** identifies incomplete bonding
evidence. References are listed below.

For neutral/cationic Xe, the selected construction convention assigns #v2 to each
Xe=O, #m<2m> to even fluoride sets, and #v to a localized Xe–F unit as in XeF+,
combined as needed. For the retained monocations, one Xe–F contributes #v and the
remaining even F set contributes #m. Source-specific localized bonds are retained.
For the selected anions, all F contribute to #m, with #v2 for each Xe=O.

| Atom descriptions | Evidence status and example | Source |
| --- | --- | --- |
| He; Ne; Ar; Kr#n4; Xe#n4; Rn#n4 | Supported: free atoms, retaining the user's He/Ne/Ar zero-lone-pair convention | HW 461 |
| Kr#n3#m2; Xe#n3#m2 | Supported: KrF2/XeF2, three-center/four-electron description; Kr mapping by analogy | HW 465, 469; C64, pp. 1445–1448 |
| Xe#n2#m4; Xe#n#m6 | Supported by XeF4/XeF6 | HW 467–468; C64, pp. 1445–1448 |
| Kr#c+#n3#v; Xe#c+#n3#v | Supported: KrF+, XeF+, XeCl+ | HW 466, 469 |
| Kr#c+#n3#v#t; Xe#c+#n3#v#t | Supported: nitrile donor adducts; HW explicitly gives HCN→KrF+ | HW 466, 469, 473; B07, p. 3608 |
| Xe#c+#n2#v#m2 | Supported: M73 explicitly describes axial three-center bonding, an equatorial electron-pair bond and two Xe lone pairs in XeF3+. K14 independently supports C6F5XeF2+ | HW 467, 473; M73, pp. 1714, 1716; K14, pp. 11656–11657 |
| Xe#c+#n#v#m4 | Supported: G10 explicitly describes XeF5+ with two equatorial three-center/four-electron bonds, one covalent axial Xe–F bond and one Xe lone pair | HW 468; G10, p. 265 |
| Xe#c+#n#v#t#m4 | Supported as a constructed donor–acceptor description: H10 describes Xe2F11+ as two XeF5+ units sharing a fluoride ion; each Xe accepts one pair from that F−, retaining the G10 core description | HW 468; G10, p. 265; H10, pp. 3504, 3513 |
| Xe#c-#n2#m5 | Supported: XeF5−, all five equivalent F in the multicenter unit, with two Xe lone pairs | HW 467; C91, pp. 3351, 3360 |
| Xe#c-#n#m7 | Inferred, selected: XeF7− with one Xe lone pair and all seven F in #m7. P64/P66 establish the salts; E96/G10 establish structures. The Wiberg comparison below supports the selection | HW 468; P64; P66; E96, pp. 1123–1124; G10, pp. 265–266 |
| Xe#c-2#n#m8 | Inferred, selected: XeF8²− with one Xe lone pair and all eight F in the multicenter description, by analogy with XeF2/XeF+ and the D4d geometry. A96 establishes the species, geometry and inactive lone pair; the calculation check below supports the lone pair | HW 468; A96, pp. 400–401 |
| Xe#n#v6; Xe#n#v6#t | Supported: XeO3 and its 1:1 nitrile adduct; B72 explicitly allows Xe=O, and G16 uses Xe=O bonds, one Xe lone pair and N-lone-pair donation | HW 469–470; B72, p. 1127; G16, pp. 13780–13782 |
| Xe#v8 | Supported by XeO4: B72 explicitly permits four Xe=O bonds; G70 favors the double-bond description from bond lengths. Perxenate [XeO6]4− is also selected, using the two-double/four-single encoding as an accepted assumption | HW 470–471; B72, p. 1127; G70, pp. 814–815 |
| Xe#n2#v2#m2 | Supported as a formal description of XeOF2: B07 includes Xe=O valence structures and explicitly describes F–Xe–F as three-center/four-electron | HW 471; B07, pp. 3601, 3607 |
| Xe#n#v2#m4 | Supported as a formal description of XeOF4: HW supplies Xe=O; W75 explicitly uses Rundle three-center/four-electron fluoride bonding in this series | HW 471; W75, p. 1477; H11, p. 666 |
| Xe#n#v4#m2 | Supported as a formal description of XeO2F2: P06 supports the two Xe=O contributors; W75 supplies the three-center/four-electron F–Xe–F description | HW 471; W75, p. 1477; P06, pp. 1529, 1532; SI p. S2 |
| Xe#v6#m2 | Inferred: XeO3F2 with three Xe=O and an F–Xe–F multicenter unit. CH71 establishes geometry and W75 supports the fluoride assignment; the Xe=O assignments use the construction convention | HW 471; CH71, pp. 1506–1507; W75, p. 1477 |
| Xe#v4#m4 | Inferred: XeO2F4 with two Xe=O and four multicenter F contributions. H82 establishes detection and volatility; the bonding assignment uses the construction convention | HW 465, 471; H82, pp. 687–688 |
| Xe#c+#n#v5 | Supported as a formal description of XeO2F+: P06 gives two Xe=O bonds, one Xe–F bond and one Xe lone pair; its bonding analysis includes this Lewis structure | HW 471; P06, pp. 1519, 1531 |
| Xe#c+#n#v3#m2 | Supported: M93 explicitly assigns axial three-center/four-electron F–Xe–F and equatorial two-center/two-electron Xe–F, and discusses the Xe=O multiple bond in XeOF3+ | HW 471; M93, pp. 388, 391 |
| Xe#c-#n2#v2#m3 | Inferred, selected: XeOF3− with Xe=O, two Xe lone pairs and all three F in #m3. B10 distributes fluoride bonding over all three F; the Wiberg comparison below supports this construction | HW 471; B10, p. 10940; SI p. S9 |
| Xe#c-#n#v2#m5 | Supported as a formal description of XeOF5−: one Xe=O, one Xe lone pair and all five equatorial F in a six-center/ten-electron unit | HW 471; C95, pp. 1873–1874; E95, pp. 1586–1587; B00, p. 62 |
| Xe#c-#n#v4#m3 | Inferred, selected: XeO2F3− with two Xe=O, one Xe lone pair and all three F in #m3. C88 establishes the cis assignment; the Wiberg comparison below supports the selection | HW 471; C88, pp. 3767–3768 |
| F#n3#m; N#v3#d | Supported: terminal F in the multicenter fluorides; nitrile donor N | HW 465–466, 469; C64; B07; G16 |
| F#c-#n2#d2 | Inferred, selected: bridging F− donates one pair to each of two XeF5+ units in Xe2F11+ | HW 468; H10, pp. 3504, 3513; PW |

#### Primary sources read

PDFs are retained in [materials/valence](../materials/valence/); text extracts and
page images are temporary scratch intermediates.

- **C64:** Coulson, “The Nature of the Bonding in Xenon Fluorides,” *J. Chem. Soc.* (1964), 1442–1454. [DOI](https://doi.org/10.1039/JR9640001442).
- **P64:** Peacock, Selig and Sheft, “Complex Fluoroxenates(VI),” *Proc. Chem. Soc.* (1964), 285. [PDF](../materials/valence/Peacock-1964-fluoroxenates.pdf).
- **P66:** Peacock, Selig and Sheft, “Complex Compounds of Xenon Hexafluoride with the Alkali Fluorides,” *J. Inorg. Nucl. Chem.* **28** (1966), 2561–2567. [DOI](https://doi.org/10.1016/0022-1902(66)80380-9).
- **G70:** Gundersen, Hedberg and Huston, “Molecular Structure of Xenon Tetroxide, XeO4,” *J. Chem. Phys.* **52** (1970), 812–815. [DOI](https://doi.org/10.1063/1.1673060).
- **CH71:** Claassen and Huston, “Matrix-Isolation Raman and Infrared Spectra, Molecular Structure, and Force Constants for XeO3F2,” *J. Chem. Phys.* **55** (1971), 1505–1507. [DOI](https://doi.org/10.1063/1.1676271).
- **B72:** Bartlett et al., “The Crystal Structure of Xenon(II) Fluoride Fluorosulfate, FXeOSO2F,” *Inorg. Chem.* **11** (1972), 1124–1127. [DOI](https://doi.org/10.1021/ic50111a042).
- **M73:** McKee, Zalkin and Bartlett, “Crystal Structure of [XeF3+][Sb2F11−],” *Inorg. Chem.* **12** (1973), 1713–1717. [DOI](https://doi.org/10.1021/ic50126a001); [PDF](../materials/valence/McKee-1973-XeF3-cation.pdf).
- **W75:** Willett, LaBonville and Ferraro, “Normal Coordinate Treatment of XeO2F2,” *J. Chem. Phys.* **63** (1975), 1474–1478. [DOI](https://doi.org/10.1063/1.431510).
- **H82:** Huston, “Chemical and Physical Properties of Some Xenon Compounds,” *Inorg. Chem.* **21** (1982), 685–688. [DOI](https://doi.org/10.1021/ic00132a043).
- **S85:** Stein, “New Evidence that Radon is a Metalloid Element: Ion-exchange Reactions of Cationic Radon,” *J. Chem. Soc., Chem. Commun.* (1985), 1631–1632. [DOI](https://doi.org/10.1039/C39850001631).
- **C88:** Christe and Wilson, “New Syntheses and Properties of XeO2F2, Cs+XeO2F3−, and NO2+[XeO2F3·nXeO2F2]−,” *Inorg. Chem.* **27** (1988), 3763–3768. [DOI](https://doi.org/10.1021/ic00294a019).
- **C91:** Christe et al., “The Pentafluoroxenate(IV) Anion, XeF5−: The First Example of a Pentagonal Planar AX5 Species,” *JACS* **113** (1991), 3351–3361. [DOI](https://doi.org/10.1021/ja00009a021).
- **M93:** Mercier et al., “The Oxotrifluoroxenon(VI) Cation: X-ray Crystal Structure of XeOF3+SbF6− and a Solution 17O and 129Xe Nuclear Magnetic Resonance Study of the 17,18O-Enriched XeOF3+ Cation,” *Inorg. Chem.* **32** (1993), 386–393. [DOI](https://doi.org/10.1021/ic00056a009).
- **C95:** Christe et al., “On the Structure of the [XeOF5]− Anion and of Heptacoordinated Complex Fluorides Containing One or Two Highly Repulsive Ligands or Sterically Active Free Valence Electron Pairs,” *Inorg. Chem.* **34** (1995), 1868–1874. [DOI](https://doi.org/10.1021/ic00111a039).
- **E95:** Ellern and Seppelt, “XeOF5−, an Anion with Pentagonal-Pyramidal Structure,” *Angew. Chem. Int. Ed. Engl.* **34** (1995), 1586–1587. [DOI](https://doi.org/10.1002/anie.199515861). Reports five approximately coplanar F and apical O; assigns partial Xe–O double-bond character and a trans lone pair.
- **A96:** Adam, Ellern and Seppelt, “Structural Principles of the Coordination Number Eight: WF8²−, ReF8²−, and XeF8²−,” *Chem. Eur. J.* **2** (1996), 398–402. [DOI](https://doi.org/10.1002/chem.19960020408); [PDF](../materials/valence/Adam-1996-XeF8-structure.pdf).
- **E96:** Ellern, Mahjoub and Seppelt, “Structures of XeF7− and Xe2F13−,” *Angew. Chem. Int. Ed. Engl.* **35** (1996), 1123–1125. [DOI](https://doi.org/10.1002/anie.199611231).
- **B00:** Baran, “Mean Amplitudes of Vibration of the Pentagonal Pyramidal XeOF5− and IOF5²− Anions,” *J. Fluorine Chem.* **101** (2000), 61–63. [DOI](https://doi.org/10.1016/S0022-1139(99)00194-3).
- **P06:** Pointner, Suontamo and Schrobilgen, “Syntheses and X-ray Crystal Structures… and Computational Studies of the XeO2F+ and FO2XeFXeO2F+ Cations and Related Species,” *Inorg. Chem.* **45** (2006), 1517–1534. [DOI](https://doi.org/10.1021/ic0510260); [SI](https://doi.org/10.1021/ic0510260.s002), p. S2 and Table S8.
- **B07:** Brock et al., “XeOF2, F2OXeN≡CCH3, and XeOF2·nHF: Rare Examples of Xe(IV) Oxide Fluorides,” *JACS* **129** (2007), 3598–3611. [DOI](https://doi.org/10.1021/ja0673480).
- **B10:** Brock et al., “XeOF3−, an Example of an AX3YE2 Valence Shell Electron Pair Repulsion Arrangement…,” *JACS* **132** (2010), 10935–10943. [DOI](https://doi.org/10.1021/ja103730c); [SI](https://doi.org/10.1021/ja103730c.s001), p. S9, Scheme S1.
- **G10:** Grant et al., “Heats of Formation of XeF3+, XeF3−, XeF5+, XeF7+, XeF7−, and XeF8 from High Level Electronic Structure Calculations,” *Inorg. Chem.* **49** (2010), 261–270. [DOI](https://doi.org/10.1021/ic901956g); [SI](../materials/valence/Grant-2010-xenon-fluoride-ions-SI.pdf), Table SM-3.
- **H10:** Hughes, Mercier and Schrobilgen, “Syntheses, Raman Spectra, and X-ray Crystal Structures of [XeF5][μ-F(OsO3F2)2] and [M][OsO3F3] (M = XeF5+, Xe2F11+),” *Inorg. Chem.* **49** (2010), 3501–3515. [DOI](https://doi.org/10.1021/ic100062x).
- **H11:** Hughes et al., “A Raman Spectroscopic Study of the XeOF4/XeF2 System and the X-ray Crystal Structure of α-XeOF4·XeF2,” *J. Fluorine Chem.* **132** (2011), 660–668. [DOI](https://doi.org/10.1016/j.jfluchem.2011.05.010).
- **P13:** Papakondylis, “Xenon oxides, sulfides, and oxysulfides,” *Comput. Theor. Chem.* **1015** (2013), 15–20. [DOI](https://doi.org/10.1016/j.comptc.2013.03.030).
- **K14:** Koppe et al., “Xenon(IV)–Carbon Bond of [C6F5XeF2]+; Structural Characterization and Bonding…,” *Inorg. Chem.* **53** (2014), 11640–11661. [DOI](https://doi.org/10.1021/ic501831j).
- **G16:** Goettel et al., “Syntheses and Structures of Xenon Trioxide Alkylnitrile Adducts,” *Angew. Chem. Int. Ed.* **55** (2016), 13780–13783. [DOI](https://doi.org/10.1002/anie.201607583).

Full text not retrieved (2026-10-31):

- **P88:** Popov, A. I.; Fadeeva, N. E.; Kiselev, Yu. M., “Oxidimetric determination of active fluorine in fluorides used as fluorinating and oxidizing agents,” *Zhurnal Analiticheskoi Khimii* **43**(3) (1988), 465–471. User-supplied reference from the XeF8²− compound searches; no online full text located, bonding content unassessed.

#### XeF8²− calculation check (2026-10-31)

PySCF 2.14.0, RHF/def2-TZVPP with the 28-electron Xe def2 ECP, at the
user-supplied MP2 geometry: D4d, Xe–F = 2.04889054 Å, upper-square polar angle
57.32284694°, lower square reflected in z and rotated 45°. XeF2 (linear) and
XeF+ references use that same distance without optimization. Wiberg indices use
the spin-summed density in PySCF's orthonormal NAO basis.

| Species | Xe–F Wiberg index per contact | Sum at Xe |
| --- | --- | --- |
| XeF8²− | 0.38480 | 3.07842 |
| XeF2 | 0.48208 | 0.96417 |
| XeF+ | 0.93872 | 0.93872 |

XeF8²−: 298 basis functions, 100 explicit electrons; converged RHF energy
−1123.50405686658 Eh. HOMO/LUMO energies agree with the supplied TURBOMOLE
output within 1 × 10⁻⁶ Eh. PySCF's Xe NAO charge is +3.70911, versus +3.66218
in TURBOMOLE; the population analyses are not numerically identical.
TURBOMOLE gives the Xe 5s NAO an occupation of 1.99841 electrons. The user selected
Xe#c-2#n#m8 by analogy with XeF2/XeF+ and the D4d geometry; the indices are supporting
comparisons, not a rule assigning integer #v/#m values.

#### Remaining-anion calculation checks (2026-10-31)

The same RHF/def2-TZVPP/ECP28 and PySCF NAO Wiberg calculation gives:

| Anion | Geometry used | Xe–F Wiberg indices | Xe–O Wiberg indices |
| --- | --- | --- | --- |
| XeF7− | G10, SI Table SM-3, C3v CCSD(T)/aVTZ | cap 0.36048; near triangle 0.43215 × 3; far triangle 0.47983 × 3 | — |
| XeOF3− | B10, Table 3, C2v MP2/aug-cc-pVTZ(-PP) | Fa 0.39454 × 2; Fb trans to O 0.20643 | 1.11069 |
| XeO2F3− | New converged MP2/def2-TZVPP geometry | axial 0.39019, 0.39016; equatorial 0.24292 | O trans to F 1.06450; other O 1.13495 |

E96 establishes a C3v crystal structure of CsXeF7. G10's gas-phase C2v and C3v
structures differ by only 0.02 kcal/mol; the calculation here samples C3v only.
B10 establishes XeOF3− by Raman spectroscopy and calculation, without a crystal
structure. C88 assigns cis XeO2F3− from vibrational spectra; no atomic coordinates
were found in the reviewed sources.

Published geometry parameters reconstructed for these checks (Å, degrees):
XeF7− cap r = 2.0441; near triangle r = 2.0046, polar angle 75.4; far triangle
r = 1.9667, polar angle 131.3, azimuth offset 60. XeOF3− r(Xe–O) = 1.821,
r(Xe–Fa) = 2.060, r(Xe–Fb) = 2.129, O–Xe–Fa = 94.8, O–Xe–Fb = 180.
The latter reconstruction gives Fa–Xe–Fa = 170.4°, versus the independently
rounded 170.3° in B10. All three RHF calculations converged; density normalization,
NAO orthogonality and closed-shell idempotency checks passed.
B10 SI Table S7 labels its quantities NBO bond orders, without identifying them
as Wiberg indices; they are not used as the numerical reference here.

XeO2F3− uses an unconstrained MP2/def2-TZVPP optimization with Xe ECP28 and all
explicit electrons correlated, starting from C88's cis topology with estimated
distances. It converged in 11 evaluations: E = −778.15641465 Eh, maximum gradient
component 0.00004959 Eh/bohr; no frequency calculation was performed. Coordinates
used (Å, rigidly aligned; O trans to F first, equatorial F last):

```text
Xe   0.00000000   0.00000000   0.00000000
O    1.76385519   0.00000000   0.00000000
O   -0.28197632   1.70408983   0.00000000
F   -0.11278732  -0.11252632   1.99680396
F   -0.11285876  -0.11151882  -1.99689469
F   -2.07556951  -0.26762219   0.00005425
```

Compared with the preliminary calculation at evaluation 7, the largest Wiberg-index
change is 0.00276. The selected registry replacements put all F contributions in
#m and retain #v2 per Xe=O.

#### Selections applied

Perxenate [XeO6]4− retains Xe#v8 as an accepted assumption. Its existing
[conformance input](../umol-graph/tests/resolution/data/inorganic_ions/xeo64-.edn)
uses two Xe=O and four Xe–O− bonds, with neutral Xe. This settles the final
noble-gas example; the atom state already has independent support from XeO4.

F#c-#n2#d2 is added for the Xe2F11+ bridge, donating one pair to each
Xe#c+#n#v#t#m4 center. Its lone-pair count is left for resolution in the
[conformance input](../umol-graph/tests/resolution/data/donor_acceptor/xe2f11+.edn).

| Implemented replacement | Replaced entry | Example |
| --- | --- | --- |
| Xe#c-#n2#m5 | Xe#c-#n2#v#m4 | XeF5− |
| Xe#c-#n#v2#m5 | Xe#c-#n#v3#m4 | XeOF5− |
| Xe#c-2#n#m8 | Xe#c-2#n#v2#m6 | XeF8²− |
| Xe#c-#n#m7 | Xe#c-#n#v#m6 | XeF7− |
| Xe#c-#n2#v2#m3 | Xe#c-#n2#v3#m2 | XeOF3− |
| Xe#c-#n#v4#m3 | Xe#c-#n#v5#m2 | XeO2F3− |

These mappings retain the formal charge on Xe and use the existing F#n3#m
description. Each n-fluoride multicenter unit contains 2n electrons. The states
conserve electrons under the registry defaults; multicenter descriptions remain
outside the Counts comparison.

Rn#c+2#n3 is excluded by the user: S85, pp. 1631–1632, proposes Rn²⁺ but does
not distinguish it from RnF+. The default entry and ionic RnF2 conformance fixture
are removed. Rn#n3#m2, Rn#n2#m4 and Rn#n#m6 are also removed for insufficient
support for their compound descriptions (HW 465, 469; S85); only Rn#n4 remains.

Xe#m8 is removed: its oxide/perxenate examples lack support for an all-multicenter
description (HW 470–471; P13). Neutral XeF8 is unestablished in HW and was not an
example for this entry. Obsolete oxide examples are removed from retained fluoride
entries; B07 supersedes HW's polymeric XeOF2 interpretation.

Xe#t#m8 and Xe#v8#t and their MeCN→XeO4 fixtures are removed by the user's
decision: no verified XeO4 adduct was found; G16, p. 13780, eq. 2, instead gives
MeCN→XeO3 + ½O2. N#v3#d is retained for supported N-donor complexes.

MeCN→XeO3 retains only Xe#n#v6#t, supported by G16. The conflicting
Xe#n#t#m6 entry is removed.

#### Further extracted examples and qualifications

| Examples | Source detail | HW pages |
| --- | --- | --- |
| Xe2F3+, Kr2F3+, Xe2F11+, Xe2F13−; oligomeric XeF6 | Bridging/association; Kr2F3+ is asymmetric | 466, 468–469 |
| XeCl2, XeCl4, XeBr2, XeClF | XeCl2 in a cold condensate; XeCl2/4 and XeBr2 as decay products; XeClF mentioned | 465, 469 |
| HeH+, EgH+, HArF, HKrCl, HKrCN, HXeCl/Br/I/CN/NC/C2H, XeH2 | Gas-phase and matrix chemistry; HW describes HArF as H–Ar+ / F− | 472 |
| Ne₂⁺, Ar₂⁺, Kr₂⁺, Xe₂⁺ | Mass spectrometry; an Xe₂⁺ salt is also reported | 472 |
| Kr(OTeF5)2; Xe(OTeF5)2/4/6; XeF(OAc), Xe(OAc)2; XeF[N(SO2F)2], Xe[N(SO2F)2]2 | Further O/N ligand examples | 472–473 |
| C6F5Xe+ and MeCN adduct; C6F5XeF/Cl/CN, C6F5XeOCO(C6F5), Xe(C6F5)2, C6F5XeF2+ | Xe–C examples | 473 |

Not selected: Xe#n3#d for AuXe4²+/Au–Xe/Hg–Xe complexes (HW 473);
Xe#c+#n3#m for Xe₂⁺ (trial encoding, HW 472); Xe#c+2#n3#t2 for
Xe(FH)2²+ (proposed interpretation in HW 466).

He#n, Ne#n4 and Ar#n4 are excluded by the user: retain the zero-lone-pair
descriptions and current valence limits. The trial descriptions He#c+#v,
Ne#c+#n3#v, Ar#c+#n3#v, Ne#n3#m2, Ar#n3#m2 and Ne#c+#n3#m are not selected.

HW marks ArF2, KrF4, neutral XeF8 and XeF7+ as unestablished; ArF+ and XeOF6
are questioned. XeO, XeO2 and Xe(OH)2 occur as possible intermediates.

#### Implementation checks

F#n3#m, F#c-#n2#d2, N#v3#d and the selected noble-gas descriptions are in the default registry.
All noble-gas entries carry HW volume/page references. After the selections above,
the noble-gas increment retains 26 conformance inputs, including
[XeO3](../umol-graph/tests/resolution/data/inorganic_small/xeo3.edn),
[XeOF4](../umol-graph/tests/resolution/data/inorganic_small/xeof4.edn),
[XeF5−](../umol-graph/tests/resolution/data/inorganic_ions/xef5-.edn),
[XeOF5−](../umol-graph/tests/resolution/data/inorganic_ions/xeof5-.edn),
[XeF8²−](../umol-graph/tests/resolution/data/inorganic_ions/xef8-2.edn) and
[XeF7−](../umol-graph/tests/resolution/data/inorganic_ions/xef7-.edn).
The existing [XeOF3−](../umol-graph/tests/resolution/data/inorganic_ions/xeof3-.edn)
and [XeO2F3−](../umol-graph/tests/resolution/data/inorganic_ions/xeo2f3-.edn) inputs
also use the selected replacements. Nine further fixtures cover KrF2, KrF+,
HCN→KrF+, XeF4, XeF6, XeF+, XeF3+, XeF5+ and Xe2F11+. All **33 retained
noble-gas states** and the three F/N companions are now exercised. Both AtomTyping
policies match the specified states, including two lone pairs on the bridging F−
and one on each Xe in Xe2F11+.
All **711 conformance cases and 366 valence unit tests pass**,
preserving the 685 cases preceding the HW increment.
The perxenate encoding resolves under both policies with AtomTyping and Counts.
No noble-gas selections remain open.

Earlier Counts probes with explicit H=0 matched the XeO3, XeO4 and XeO6⁴−
encodings. Multicenter/dative cases remain outside its supported scope.
ValenceTable is unchanged.

### Halogens

Source: [HW 103, Band 1, chapter XII](<../materials/valence/Holleman-Wiberg - Anorganische Chemie - 103. Auflage - Band 1.pdf#page=513>),
printed pp. 475–552, PDF pp. 513–590. References below are printed HW I pages.
Descriptions use AtomTypeRegistry::raise_defaults; the formulas are source examples,
and the atom strings are constructed encodings. The user's L-table selections are
recorded below; **implementation is paused for clarification of overlapping descriptions**.
The M and D tables remain unreviewed except for the F entry in D04, now implemented
for Xe2F11+ in the noble-gas slice.

#### Existing coverage

The registry has 37 halogen descriptions: F 6, Cl 11, Br 9, I 9, At 2.

| Existing descriptions | Examples / source qualification | HW pages |
| --- | --- | --- |
| X#n3#v; X#n3#u (X = F, Cl, Br, I, At); X#n3#h (X = F, Cl, Br, I) | X2, HX, interhalogens, terminal halogens and free atoms; At chemistry rests on tracer experiments | 477, 483, 486, 489–490, 496–506 |
| X#c-#n4 (X = F, Cl, Br, I) | Halide salts | 493, 498–505 |
| X#c+#n3; X#c+#n2#v2 (X = Cl, Br, I) | Atomic cations; X3+ and XF2+ salts. HW distinguishes free gas-phase I+ from condensed-phase complexes | 491–492, 510 |
| X#n2#v3; X#n#v5; X#v7 (X = Cl, Br, I) | XF3/XF5, IF7; oxoacids/oxoanions, oxides and fluoride oxides; I(OTeF5)3/5 | 508–510, 512–546 |
| Cl#n#u#v4; Cl#u#v6 | ClO2; transient ClO3 | 531, 536–540 |
| F#n3#m | Multicenter terminal F, added in the noble-gas slice | 465, 469 |
| F#c-#n2#d2 | Bridging F− in Xe2F11+, added in the noble-gas slice; covers the F entry in D04 | 468 |

#### Descriptions and selection status

Each X expands only over the elements specified in its row. Examples in these tables
also supply candidates for conformance inputs. The user accepted L01–L09 and L11,
with separate At oxoanion entries required in L03. L10 is covered by the user's
N#c-2#n3#v addition. These selections do not approve adding both representations of
the same examples. The following overlaps were not resolved in the extraction:

| Rows requiring clarification | Overlap / source uncertainty |
| --- | --- |
| L04 fluoride examples | The cited ion list establishes XF4+/XF6+ species; the extraction did not verify their localized bond description. Oxide cations and I(OH)6+ are separate examples in this row |
| L05 ↔ M02 | Alternative descriptions of X3− and XY2−. HW explicitly discusses 3-center/4-electron bonding for I3− on p. 494; L05 was not independently justified by that example |
| L06 ↔ M03, M04, M11, M12 | Alternative descriptions of the same polyhalides and oxofluorides. H2IO4− is a separate example for I#c-#n#v6 |
| L07 ↔ M05, M13 | IF8− overlaps M05; IOF6−/IO2F4− overlap M13. The extraction did not establish which bond description the cited source supports |
| M06–M10 ↔ existing localized entries | Alternative encodings of XF3/XF5/IF7 and fluoride oxides already represented in the registry. The specific multicenter partitions were constructed, not verified from the cited passages |

The overlapping descriptions are **unresolved alternatives**, not independent additions.

| ID | Localized descriptions | Examples / source qualification | HW pages |
| --- | --- | --- | --- |
| L01 | F#c+#n2#v2; F#c+#n2#v#h; F#c+#n2#h2 | H2F+ in liquid HF; explicit, mixed and implicit H | 498–499 |
| L02 | At#c-#n4; At#c+#n3 | Astatides; At+ precipitation as At2S and cathodic deposition; tracer chemistry | 490, 492, 499, 502 |
| L03 | At#n2#v3; At#n#v5; At#v7 | HAtO2/HAtO3/HAtO4 in the general oxoacid series; AtO3−/AtO4− in the redox diagram | 512–513 |
| L04 | X#c+#n#v4; X#c+#v6 (X = Cl, Br, I) | XF4+/XF6+ salts; ClO2+, IO2+, ClO2F2+, I(OH)6+ | 510, 528, 540, 543, 545 |
| L05 | X#c-#n3#v2 (X = Cl, Br, I, At) | X3− (Cl, Br, I); ClF2−/BrF2−/IF2−, BrCl2−, ICl2−, IBr2−, AtBr2−/AtI2− | 493–494, 510 |
| L06 | X#c-#n2#v4; X#c-#n#v6 (X = Cl, Br, I) | XF4−, BrCl4−/ICl4−/IBr4−; BrF6−/IF6−, ClO2F2−/BrO2F2−, BrOF4−, H2IO4− | 510, 527, 545 |
| L07 | I#c-#v8 | IF8−, IOF6−, cis-/trans-IO2F4− | 510–511, 546 |
| L08 | X#n2#u#v2 (X = Cl, Br, I) | XO radicals with X=O and the unpaired electron on X; transient/matrix species, HW discusses delocalization | 530–531 |
| L09 | X#n#u#v4; X#u#v6 (X = Br, I) | XO2/XO3 radicals; transient/matrix species | 531 |
| L10 | N#c-2#n3#v | N in [ClO3N]2−, retaining three Cl=O and one Cl–N; HW reports K2[ClO3N] | 546 |
| L11 | At#n3#h | HAt mentioned in the HX redox series; no preparation or structural evidence supplied there | 502 |

| ID | Multicenter descriptions | Examples / source qualification | HW pages |
| --- | --- | --- | --- |
| M01 | X#n3#m (X = Cl, Br, I) | Terminal X in X3− and mixed polyhalides; terminal F already covered | 493–494, 510–511 |
| M02 | X#c-#n3#m2 (X = F, Cl, Br, I, At) | Central X in X3− (F, Cl, Br, I), XY2−; F3− only in a 15 K matrix, decomposes at 40 K; AtBr2−/AtI2− | 493–494, 510 |
| M03 | X#c-#n2#m4 (X = Cl, Br, I) | XF4− and BrCl4−/ICl4−/IBr4− | 510–511 |
| M04 | X#c-#n#m6 (X = Br, I) | BrF6−/IF6− | 510–511 |
| M05 | I#c-#m8 | IF8− | 510–511 |
| M06 | X#n2#v#m2; X#n#v#m4 (X = Cl, Br, I) | XF3/XF5; one localized X–F, remaining F in multicenter bonds | 508–509 |
| M07 | I#v#m6 | IF7 with one localized I–F | 508–509 |
| M08 | X#n#v3#m2 (X = Cl, Br, I) | XOF3 with X=O and one localized X–F | 544–545 |
| M09 | X#v5#m2 (X = Cl, I) | ClO2F3 and monomeric IO2F3 above 100 °C, with two X=O and one localized X–F; IO2F3 forms an O-bridged dimer below 100 °C | 544–546 |
| M10 | I#v3#m4 | IOF5 with I=O and one localized I–F | 544, 546 |
| M11 | X#c-#n#v4#m2 (X = Cl, Br) | ClO2F2−/BrO2F2− with two X=O | 545 |
| M12 | Br#c-#n#v2#m4 | BrOF4− with Br=O | 545 |
| M13 | I#c-#v2#m6; I#c-#v4#m4 | IOF6−; cis-/trans-IO2F4−, with I=O | 546 |

| ID | Donor–acceptor descriptions | Examples / source qualification | HW pages |
| --- | --- | --- | --- |
| D01 | X#c+#n3#t2 (X = Br, I, At) | [X(py)2]+, two N→X bonds; py = pyridine | 491–492 |
| D02 | I#n3#v#t | Me3N→I2; HW describes a charge-transfer complex | 494–495 |
| D03 | X#c-#n3#d (X = F, Cl, Br, I, At) | Monodentate halogenido ligands, using X− donation | 499 |
| D04 | X#c-#n2#d2 (X = F, Cl, Br, I) | Bridging halogenido ligands, two donations; Nb fluoride, Al halides, Co–I–Co examples | 495–496, 499 |
| D05 | X#c-#n#d3 (X = F, Cl) | Three donations: [F(XeOF4)3]−, [F(SbF5)3]−, bridging Cl in TeCl4 | 499 |
| D06 | X#n2#v#d (X = F, Cl, Br, I) | Covalent bond plus donation in halide bridges; I2 donation in polymeric Ag(I2)+ | 489, 495–496, 499 |
| D07 | X#n2#v3#t (X = Br, I) | BrF3 polymer contacts; I acceptor in I2Cl6 with Cl donor D06 | 508–509 |
| D08 | I#n2#v3#t2; I#n#v5#t | IF3 polymer contacts; Me2SO adduct of HIO3 | 509, 527 |
| D09 | Cl#n#u#v4#t | Stable ClO2·py adduct, if encoded with N→Cl; HW does not specify the attachment in this passage | 538 |
| D10 | N#v2#a#d | Pyridine donor for D01/D09; N#v3#d for Me3N already exists | 491–492, 538 |
| D11 | O#n#v2#d | Me2SO donor for D08, using S=O and O→I | 527 |

#### Further examples and source qualifications

| Examples | Extraction result / unresolved encoding | HW pages |
| --- | --- | --- |
| Br2+/I2+, Cl4+, I4²+, Br5+/I5+, I15³+ | Isolated salts; Cl4+ decomposes above −78 °C. HW describes bond order 1.5 for X2+ and multicenter links in larger cations. Atom descriptions remain open | 490–493 |
| X2− (X = F, Cl, Br, I) | Irradiated glasses/matrices; M01/M02 (existing F#n3#m) provide a 2-center/3-electron encoding, with the unpaired electron on the bond | 493–494 |
| Br5−; higher polyiodides | Aggregates with linear, branched or extended structures; complete encodings remain open | 493–494 |
| HF2−, [HnFn+1]− (n = 1, 2, 3, 4, 5, 7) | Salts and hydrogen-bonded aggregates; encoding remains open | 498 |
| Br2·toluene; F2/Cl2/Br2/ClF/BrCl with H2O, H2S or NH3 | Charge-transfer complex or short-lived prereactive associates; no dative description assigned solely from the contact | 494–495 |
| ClO4, Cl2O2+; polymeric IO+, IO2+, I3O6+ and iodine oxides | ClO4 transient; Cl2O2+ salt has a π*–π* association. Local oxide units reuse existing states/L04; full association encodings remain open | 531, 539, 542–544 |
| ClF7, BrF7, F3+, I(OAc)7; I7+, I(OSO2F)5, XO3+ | First four unestablished in HW; last three possible or questioned | 491–492, 508–510, 545 |

#### Registry and ValenceTable comparison

The initial extraction contained **92 distinct candidate descriptions** (89 halogen,
3 N/O companions), none then present under the registry's loading defaults; L10 has
since been added by the user. All pass electron conservation. This check does not
establish source support for a bond description or settle the alternatives above.
All 34 localized descriptions pass Counts classification with their literal H counts;
no ValenceTable exclusion was found for this set, including I#c-#v8. Multicenter and
dative descriptions are outside this Counts comparison.

Counts uses saturation targets, not an allowed-valence whitelist. Both policies
resolve H2F+, ClF4+, localized I3−/IF8− and BrO with H=0 in probes. With both H and
spin unspecified, BrO is underdetermined under Strict and receives H under
MostSaturated; supplying the radical or H=0 retains BrO. No registry, ValenceTable
or conformance changes have been made for this chapter.
