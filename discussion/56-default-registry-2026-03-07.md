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
| H | yes | yes (+1) | yes | yes (`atoms/h.toml`, `atomic_ions/h+1.toml`, hydrides) | done |  |
| He | yes | no | n/a | yes (`atoms/he.toml`) | done | |
| Li | yes | yes (+1) | no | yes (`atoms/li.toml`, `atomic_ions/li+1.toml`) | done |  |
| Be | yes | yes (+2) | no | yes (`atoms/be.toml`, `atomic_ions/be+2.toml`) | done |  |
| B | yes | no | no | yes (`atoms/b.toml`) | done |  |
| C | yes | no | no | yes (`atoms/c.toml`) | done |  |
| N | yes | yes (-3) | no | yes (`atoms/n.toml`, `atomic_ions/n-3.toml`) | done |  |
| O | yes | yes (-2) | no | yes (`atoms/o.toml`, `atomic_ions/o-2.toml`) | done |  |
| F | yes | yes (-1) | no | yes (`atoms/f.toml`, `atomic_ions/f-1.toml`, hydrides) | done |  |
| Ne | yes | no | n/a | yes (`atoms/ne.toml`) | done |  |
| Na | yes | yes (+1) | no | yes (`atoms/na.toml`, `atomic_ions/na+1.toml`) | done |  |
| Mg | yes | yes (+2) | no | yes (`atoms/mg.toml`, `atomic_ions/mg+2.toml`) | done |  |
| Al | yes | yes (+3) | no | yes (`atoms/al.toml`, `atomic_ions/al+3.toml`) | done |  |
| Si | yes | no | no | yes (`atoms/si.toml`) | done |  |
| P | yes | yes (-3) | no | yes (`atoms/p.toml`, `atomic_ions/p-3.toml`) | done |  |
| S | yes | yes (-2) | no | yes (`atoms/s.toml`, `atomic_ions/s-2.toml`) | done |  |
| Cl | yes | yes (-1) | no | yes (`atoms/cl.toml`, `atomic_ions/cl-1.toml`) | done |  |
| Ar | yes | no | n/a | yes (`atoms/ar.toml`) | done |  |
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
| Kr | yes | no | n/a | yes (`atoms/kr.toml`) | done |  |
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
| Xe | yes | no | partial | yes (`atoms/xe.toml`) | done |  |
| Cs | yes | yes (+1) | no | yes (`atoms/cs.toml`, `atomic_ions/cs+1.toml`) | done |  |
| Ba | yes | yes (+2) | no | yes (`atoms/ba.toml`, `atomic_ions/ba+2.toml`) | done |  |
| Tl | yes | yes (+1,+3) | no | yes (`atoms/tl.toml`, `atomic_ions/tl+1.toml`, `atomic_ions/tl+3.toml`) | done |  |
| Pb | yes | yes (+2,+4) | no | yes (`atoms/pb.toml`, `atomic_ions/pb+2.toml`, `atomic_ions/pb+4.toml`) | done |  |
| Bi | yes | yes (+3) | no | yes (`atoms/bi.toml`, `atomic_ions/bi+3.toml`) | done |  |
| Po | yes | yes (-2) | no | yes (`atoms/po.toml`, `atomic_ions/po-2.toml`) | done |  |
| At | yes | yes (-1) | no | yes (`atoms/at.toml`, `atomic_ions/at-1.toml`) | done |  |
| Rn | yes | no | no | yes (`atoms/rn.toml`) | done |  |
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
| Carbonyls | Aldehydes, ketones | — | todo | |
| Carboxylic acids / derivatives | Acids, esters, amides | — | todo | |
| Organic halides | R–F, R–Cl, R–Br, R–I | — | todo | |
| Other | Amines, thiols, etc. | — | todo | |

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

Use resolution data categories:

- `atoms/` — free atoms (ground and valence/excited states).
- `atomic_ions/` — single-atom ions (discrete charge states).
- `hydrides/` — one-atom graphs with implicit H (`?{EHn}`, `implicit_h=true`); no explicit H atoms.
- `inorganic_small/` — small neutral inorganic molecules (including sextet species); explicit H only (`implicit_h=false`). Split from former `covalent/`.
- `inorganic_ions/` — small molecular ions (e.g. NH4+, OH3+); explicit H only. Split from former `covalent/`.
- `hydrocarbons/` — C/H only, implicit H (`?{CHn}`, `implicit_h=true`).
- `functional_groups/` — organic molecules with heteroatom functional groups; implicit H (`?{CHn}`, `?{OH}`, `?{NH2}`, `?{SH}`, etc., `implicit_h=true`).

Compromise: hydrides use implicit H (one-atom graph); inorganic_small, inorganic_ions, and hydrocarbons use explicit H where present, except hydrocarbons where H is in the query.

## Current Conformance Set

- `atoms/`: `h`, `he`, `li`, `be`, `b`, `c`, `n`, `o`, `f`, `ne`, `na`, `mg`, `al`, `si`, `p`, `s`, `cl`, `ar`, `k`, `ca`, `sc`, `ti`, `v`, `cr`, `mn`, `fe`, `co`, `ni`, `cu`, `zn`, `ga`, `ge`, `as`, `se`, `br`, `kr`, `rb`, `sr`, `y`, `zr`, `nb`, `mo`, `tc`, `ru`, `rh`, `pd`, `ag`, `cd`, `in`, `sn`, `sb`, `te`, `i`, `xe`, `cs`, `ba`, `tl`, `pb`, `bi`, `po`, `at`, `rn`, `la`, `ce`, `pr`, `nd`, `pm`, `sm`, `eu`, `gd`, `tb`, `dy`, `ho`, `er`, `tm`, `yb`, `lu`, `hf`, `ta`, `w`, `re`, `os`, `ir`, `pt`, `au`, `hg`, `fr`, `ra`, `ac`, `th`, `pa`, `u`, `np`, `pu`, `am`, `cm`, `bk`, `cf`, `es`, `fm`, `md`, `no`, `lr`, `rf`, `db`, `sg`, `bh`, `hs`
- `atomic_ions/`: `h+1`, `li+1`, `be+2`, `mg+2`, `al+3`, `n-3`, `o-2`, `f-1`, `na+1`, `p-3`, `s-2`, `cl-1`, `k+1`, `ca+2`, `sc+3`, `ti+3`, `ti+4`, `v+2`, `v+3`, `cr+2`, `cr+3`, `mn+2`, `mn+3`, `fe+2`, `fe+3`, `co+2`, `co+3`, `ni+2`, `cu+1`, `cu+2`, `zn+2`, `ga+3`, `as-3`, `se-2`, `br-1`, `rb+1`, `sr+2`, `y+3`, `zr+4`, `nb+3`, `mo+2`, `mo+3`, `ru+2`, `ru+3`, `rh+2`, `rh+3`, `pd+2`, `ag+1`, `ag+2`, `cd+2`, `in+3`, `sb-3`, `te-2`, `i-1`, `cs+1`, `ba+2`, `tl+1`, `tl+3`, `pb+2`, `pb+4`, `bi+3`, `po-2`, `at-1`, `la+3`, `ce+3`, `ce+4`, `pr+3`, `pr+4`, `nd+3`, `pm+3`, `sm+2`, `sm+3`, `eu+2`, `eu+3`, `gd+3`, `tb+3`, `tb+4`, `dy+3`, `ho+3`, `er+3`, `tm+3`, `yb+2`, `yb+3`, `lu+3`, `hf+4`, `pt+2`, `au+1`, `au+3`, `hg+2`, `hg+2_dimer`, `fr+1`, `ra+2`, `ac+3`, `th+3`, `th+4`, `pa+4`, `pa+5`, `u+3`, `u+4`, `u+5`, `u+6`, `np+3`, `np+4`, `np+5`, `np+6`, `pu+3`, `pu+4`, `pu+5`, `pu+6`, `pu+7`, `am+2`, `am+3`, `am+4`, `cm+3`, `cm+4`, `bk+3`, `bk+4`, `cf+2`, `cf+3`, `cf+4`, `es+2`, `es+3`, `fm+2`, `fm+3`, `md+2`, `md+3`, `no+2`, `lr+3`
- `hydrides/`: `h2`, `bh3`, `ch4`, `nh3`, `h2o`, `hf`, `sih4`, `ph3`, `h2s`, `hcl`, `ash3`, `h2se`, `hbr`, `sbh3`, `h2te`, `hi`
- `inorganic_small/`: neutral small inorganics (explicit H only): e.g. `bf2`, `bf3`, `bn`, `c2`, `ch3`, `ch2_singlet`, `ch2_triplet`, `sih3`, `sih2_*`, `nh2`, `ph2`, `ash2`, `sbh2`, `oh`, `sh`, `seh`, `teh`, `br2`, `cl2`, `cn-`, `f2`, `i2`, `n2`, `o2`, etc. Added nitrogen oxides (NO, NO2, N2O, N2O3, N2O4), carbon compounds (CO2, CS2, COS, HCN), oxyacids/hydroxides (HNO2, HNO3, HClO, H2CO3, H3BO3, HOBr, HOI), N-H compounds (N2H4, NH2OH, N2H2), peroxides (H2O2, H2S2), halides (NF3, NCl3, PF3, PCl3, SF2, SiF4, SiCl4, BCl3, ClF, BrF), and other (O3, HNCO, HOCN, HSCN, HNCS).
- `inorganic_ions/`: molecular ions (explicit H only): e.g. `bcl4-`, `cf3-`, `ccl3+`, `nh2-`, `nh4+`, `ph2-`, `ph4+`, `ash2-`, `ash4+`, `sbh2-`, `sbh4+`, `oh-`, `oh3+`, `sh-`, `sh3+`, `seh-`, `seh3+`, `teh-`, `teh3+`, `f-`, `cl-`, `br-`, `i-`, `sif3-`, `sicl3+`, etc.
- `hydrocarbons/`: `ethane`, `propane`, `butane`, `isobutane`, `pentane`, `isopentane`, `neopentane`, `hexane`, `ethene`, `propene`, `but-1-ene`, `but-2-ene`, `2-methylpropene`, `2-methylbut-2-ene`, `ethyne`, `propyne`, `but-1-yne`, `but-2-yne`, `allene`, `buta-1,3-diene`, `penta-1,4-diene`, `butatriene`, `cyclopropane`, `cyclobutane`, `cyclopentane`, `cyclohexane`, `cycloheptane`, `cyclopropene`, `cyclobutene`, `cyclopentene`, `cyclopentadiene`, `cyclohexene`, `cyclohexa-1,3-diene`, `cyclohexa-1,4-diene`, `bicyclo-1.1.0-butane`, `bicyclo-2.1.0-pentane`, `bicyclo-2.2.0-hexane`, `bicyclo-3.3.0-octane`, `bicyclo-4.3.0-nonane`, `decalin`, `spiropentane`, `spirohexane`, `spiroheptane`, `spirononane`, `bicyclo-1.1.1-pentane`, `norbornane`, `bicyclo-2.2.2-octane`, `tetrahedrane`, `prismane`, `cubane`, `adamantane`
- `functional_groups/`: methyl/ethyl × {fluoride, chloride, bromide, iodide, azide, isocyanide, cyanate, thiocyanate, isothiocyanate, nitrite, nitrate, thiol, selenol} (26); {propyl, isopropyl, butyl, isobutyl, tert-butyl, cyclopropyl, cyclohexyl}-chloride (7); dichloromethane, chloroform, carbon-tetrachloride, 1,1-dichloroethane, 1,2-dichloroethane, 1,2-dichlorocyclohexane, 1-chloro-2-hydroxyethane (7); methanol, ethanol, propan-1-ol, propan-2-ol, butan-1-ol, butan-2-ol, 2,2-dimethylpropan-2-ol, hydroxycyclopropane, hydroxycyclohexane (9); dimethylether, diethylether, 1,2-dimethoxyethane (3); methylamine, ethylamine, dimethylamine, trimethylamine, tetramethylammonium, trimethylamine-oxide (6); nitrosomethane, nitrosoethane, nitromethane, nitroethane, nitrosocyclohexane, nitrocyclohexane (6). Total: 64 files, `implicit_h=true`.
- `hypervalent/`: expanded-octet species (all `implicit_h=false`): SO2, SO3, H2SO3, H2SO4, SF4, SF6, H3PO3, H3PO4, P4O10, PF5, PCl5, POCl3, HClO2, HClO3, HClO4, ClF3, IF5. Total: 17 files.

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